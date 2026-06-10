# Module with helper functions for xspec spectra
import os
import sys
import numpy as np
import xspec
from wurlitzer import pipes

class SpectrumLoader:
    '''This is a helper class to store the file information for a single spectrum'''
    def __init__(self, pi, path=None, back=None, rmf=None, arf=None):
        '''
        Initialize the SpectrumLoader to manage and validate files to be loaded into pyXspec.

        :param pi: The filename(s) or full path(s) of the primary spectrum file(s).
        :type pi: str or sequence of str
        :param path: Optional common directory path where the spectral files are located. 
            This path will be prefixed to all spectra and associated bundle of files (rmf/arf/back) 
            if it is a single string. If this path is a list, 
            each path will be prefixed to each spectrum's bundle of files. 
        :type path: str or sequence of str, optional
        :param back: The filename(s) or full path(s) of the background spectrum file(s).
        :type back: str or sequence of str, optional
        :param rmf: The filename(s) or full path(s) of the Redistribution Matrix File(s) (RMF).
        :type rmf: str or sequence of str, optional
        :param arf: The filename(s) or full path(s) of the Ancillary Response File(s) (ARF).
        :type arf: str or sequence of str, optional

        :raises TypeError: If any parameter does not match its expected variable type, or 
                           if an ARF is supplied for an index where the RMF is None.
        :raises FileNotFoundError: If the specified base directory or any resolved 
                                   spectral file does not exist on disk.
        :raises ValueError: If the `pi` parameter is provided as an empty sequence.
        
        ## Examples:

            ```python
            # Example 1: Loading a single spectrum with a common directory path
            loader = SpectrumLoader(path="InputData/ObsID14571/",
                pi="src_spectrum.pi",
                back="bg_spectrum.pi",
                rmf="src.rmf",
                arf="src.arf"
            )
            # To load into xspec, use:
            loader.xspec_load()


            # Example 2: Broadcasting single path across multiple PI spectra
            paths = ["Obs1", "Obs2"]
            pi_files = ["src_grp_1.pi", "src_grp_2.pi""]
            bg_files = ["bg_1.pi", "bg_2.pi"]
            rmf_files = ["rmf_1.rmf", "rmf_2.rmf"]
            loader = SpectrumLoader(
                path=paths,
                pi=pi_files,
                back=bg_files,
                rmf=rmf_files,  
                arf=None  
                )
            ```
        
        '''
        # Validate the input type
        self._validate_input_type(pi, "pi", allow_none=False)
        self._validate_input_type(path, "path", allow_none=True)
        self._validate_input_type(back, "back", allow_none=True)
        self._validate_input_type(rmf, "rmf", allow_none=True)
        self._validate_input_type(arf, "arf", allow_none=True)

        # 1. Store the base path and validate it if it's a single string
        # (a bit redundant, cause it would also fail cleanly at the
        # broadcast call below, but easier for the user if this fails first)
        self.path = path
        if isinstance(self.path, str) and not os.path.exists(self.path):
            raise FileNotFoundError(f"The specified common path directory does not exist: '{self.path}'")   

        # 2. Normalize 'pi' into a list because it defines how many spectra we are loading
        self._pi_list = [pi] if isinstance(pi, str) else list(pi)
        num_spectra = len(self._pi_list)  

        if num_spectra == 0:
            raise ValueError("You must provide at least one spectrum ('pi') file.")

        self.num_spectra = num_spectra

        # 3. Standardize all other inputs to match the length of pi_list
        self._path_list = self._broadcast_input(path, num_spectra, "path")
        self._back_list = self._broadcast_input(back, num_spectra, "back")
        self._rmf_list  = self._broadcast_input(rmf, num_spectra, "rmf")
        self._arf_list  = self._broadcast_input(arf, num_spectra, "arf")

        # 4. Process, validate, and build the clean internal lists
        self.pi = []
        self.back = []
        self.rmf = []
        self.arf = []

        for i in range(num_spectra):
            # Temporarily stage the specific path for the _resolve_and_validate method
            current_path = self._path_list[i]
            
            # Validate the path directory for this specific spectrum if it's a unique string
            if current_path and not os.path.exists(current_path):
                raise FileNotFoundError(f"The specified common path directory does not exist: '{current_path}'")

            # Resolve individual files using your strict structural logic
            resolved_pi = self._resolve_and_validate(self._pi_list[i], current_path, "Mandatory spectrum file")
            self.pi.append(resolved_pi)
            
            resolved_back = self._resolve_and_validate(self._back_list[i], current_path, "Background file") if self._back_list[i] else None
            self.back.append(resolved_back)
            
            resolved_rmf = self._resolve_and_validate(self._rmf_list[i], current_path, "RMF file") if self._rmf_list[i] else None
            self.rmf.append(resolved_rmf)
            
            resolved_arf = self._resolve_and_validate(self._arf_list[i], current_path, "ARF file") if self._arf_list[i] else None
            self.arf.append(resolved_arf)

        for i in range(num_spectra):
            # check that no arf is defined if rmf is None
            if self.rmf[i] is None:
                if self.arf[i] is not None:
                    raise TypeError(f"An arf cannot be loaded for index {i} if rmf for index {i} is None")
    
    def xspec_load(self, verbose=False):
        '''
        Load the validated spectra and their associated response files into PyXspec.

        This function wraps the low-level C++ stream execution of `xspec.Spectrum` using a system pipe 
        capture wrapper. 
        For example, pyXspec always attempts to first load the rmf/arf/back from the spectrum file's header keywords, 
        which will not work if the files are in a different directory (and pollutes the output stream with error messages), 
        and then will load the user-provided path/file. 
        Furthermore, if PyXspec encounters a critical error loading a file resquested by the user
        (e.g. if a user tried to load a spectrum file as an arf file), python will not crash.

        This function suppresses the low-level messaging, and raises clean Python exceptions when a real problem is encountered.

        :param verbose: If True, prints the underlying raw standard output and error logs 
                        received from the interactive XSpec engine. Defaults to False.
        :type verbose: bool

        :return: A single loaded XSpec Spectrum object if one file was processed, or a list 
                 of XSpec Spectrum objects if multiple files were processed.
        :rtype: xspec.Spectrum or list[xspec.Spectrum]

        :raises RuntimeError: If PyXspec fails to load a spectrum file, embedding the 
                              underlying C++ terminal crash logs.
        :raises ValueError: If a successfully loaded XSpec spectrum does not match the loader's parameters.

        ## Examples:

            ```python
            # Example 1: Loading a single spectrum with a common directory path
            loader = SpectrumLoader(path="InputData/ObsID14571/",
                pi="src_spectrum.pi",
                back="bg_spectrum.pi",
                rmf="src.rmf",
                arf="src.arf"
            )
            # To load into xspec, use:
            loader.xspec_load()


            # Example 2: Broadcasting single path across multiple PI spectra
            paths = ["Obs1", "Obs2"]
            pi_files = ["src_grp_1.pi", "src_grp_2.pi""]
            bg_files = ["bg_1.pi", "bg_2.pi"]
            rmf_files = ["rmf_1.rmf", "rmf_2.rmf"]
            loader = SpectrumLoader(
                path=paths,
                pi=pi_files,
                back=bg_files,
                rmf=rmf_files,  
                arf=None  
                )
            ```
        '''

        loaded_spectra = []

        # The content of SpectrumLoader Object has already been validated at initialization
        # so no need to do it again here. 

        # Initialize variables outside the block to catch state
        xspec_error_triggered = False
        captured_exception = None
        failed_index = None

        # 1. Start the pipe capture        
        with pipes() as (stdout, stderr):

            for i in range(len(self.pi)):
                try:
                        spec_obj = xspec.Spectrum(
                            self.pi[i],
                            backFile=self.back[i],
                            respFile=self.rmf[i],
                            arfFile=self.arf[i]
                        )
                        loaded_spectra.append(spec_obj)
                except Exception as e: # This will be triggered if the pi file exist but is of the wrong type
                    # If Python crashes here, read the 
                    xspec_error_triggered = True
                    captured_exception = e
                    failed_index = i
                    break

        if xspec_error_triggered:

            # Read whatever text is currently sitting in the capture buffer
            # Using a non-blocking read or checking the string buffer directly
            # If using standard wurlitzer, stdout.read() can hang, so we grab the buffer value:
            xspec_output = stdout.read() #if hasattr(stdout, 'read') else ""
            xspec_error = stderr.read() #if hasattr(stderr, 'read') else ""
            combined_output = f"{xspec_output}\n{xspec_error}".strip()

            if combined_output:
                if verbose:
                    print("\n" + "=" * 50)
                    print(" --- [Wrap_xspec]: XSPEC STACK LOG BEFORE CRASH --- ".center(50, "!"))
                    print("=" * 50)
                    print(combined_output)
                    print("=" * 50 + "\n")
            raise RuntimeError(
                f"[wrap_xspec]: PyXspec failed to load spectrum at index {failed_index} ('{self.pi[failed_index]}'). "
                f"Please make sure that {self.pi[failed_index]} is a valid spectrum file. "
                f"Underlying error: {captured_exception}"
                    )
        
        self._check_alldata_status(loaded_spectra)

        # 2. Normal execution path (No exceptions thrown)
        xspec_output = stdout.read()

        if verbose and xspec_output.strip():
            print("\n" + "=" * 50)
            print(" --- [Wrap_xspec]: UNDERLYING XSPEC LOADING OUTPUT --- ".center(50, "="))
            print("=" * 50)
            print(xspec_output)
            print("=" * 50 + "\n")

        print
        print('[Wrap_xspec]: Current Status of xspec.AllData:')
        print()
        xspec.AllData.show()

        # If there's only one spectrum, unpack it and return it directly.
        if len(loaded_spectra) == 1:
            return loaded_spectra[0]
            
        return loaded_spectra
 
    @staticmethod    
    def _validate_input_type(user_input, param_name, allow_none=True):
        """
        Validate that a user-supplied parameter is of an acceptable data type structure.

        Ensures that the input is either a standalone string, a valid sequence configuration 
        (list, tuple, or NumPy array) containing only strings, or None (if explicitly permitted).

        :param user_input: The raw input variable passed by the user to be type-checked.
        :type user_input: Any
        :param param_name: The name of the parameter being validated (used for error reporting).
        :type param_name: str
        :param allow_none: Whether `None` should be accepted as a valid type. Defaults to True.
        :type allow_none: bool

        :raises TypeError: If the input structure or any item within a sequence does not 
                           conform to a string or allowed types.
        
        """
        if user_input is None:
            if allow_none:
                return
            raise TypeError(f"Parameter '{param_name}' cannot be None.")

        # If it's a single item, it must be a string
        if isinstance(user_input, (str, np.str_)):
            return

        # If it's a sequence (list, tuple, numpy array), check its contents
        if isinstance(user_input, (list, tuple, np.ndarray)):
            # Catch empty sequences early
            if len(user_input) == 0:
                return
                
            for element in user_input:
                # Elements can be strings, or None (if allowed like in path lists)
                if element is None and allow_none:
                    continue
                if not isinstance(element, (str, np.str_)):
                    raise TypeError(
                        f"All elements in parameter '{param_name}' must be strings. "
                        f"Found type '{type(element).__name__}' instead."
                    )
            return

        # Explicitly reject dictionaries, sets, integers, etc.
        raise TypeError(
            f"Parameter '{param_name}' must be a string or a sequence of strings (list, tuple, numpy array). "
            f"Got type '{type(user_input).__name__}'."
        )   

    @staticmethod     
    def _to_list(user_input):
        """
        Normalize a flexible input variable into a standardized Python list.

        Converts standalone strings, None types, tuples, or NumPy arrays into a uniform, 
        iterable list format to simplify subsequent sequence processing and indexing loops.

        :param user_input: The input string, sequence, or None value to transform.
        :type user_input: Any

        :return: A standardized flat list containing the original input element(s).
        :rtype: list
        
        """
        if user_input is None or isinstance(user_input, (str, np.str_)):
            return [user_input]
        return list(user_input)

    @staticmethod     
    def _broadcast_input(user_input, expected_length, param_name):
        """
        Broadcast a scalar or sequence input to match the target spectrum dimensions.

        If the input is a scalar (a single string or None), it is duplicated into a list 
        matching `expected_length`. If the input is already an explicit sequence, this 
        method confirms that its element count strictly matches the required dimensions.

        :param user_input: The parameter value or sequence array to be broadcasted.
        :type user_input: Any
        :param expected_length: The required length of the final list (defined by the number of PI files).
        :type expected_length: int
        :param param_name: The name of the parameter being checked (used for error reporting).
        :type param_name: str

        :return: A list of length `expected_length` containing the broadcasted parameters.
        :rtype: list

        :raises ValueError: If the user explicitly provided a sequence whose length does 
                             not match the expected spectrum count.
        """

        # Scenario A: User passed None or a single string -> repeat it for every spectrum
        if user_input is None or isinstance(user_input, (str, np.str_)):
            return [user_input] * expected_length
        
        # Scenario B: User passed an explicit sequence (list, tuple, array)
        input_list = SpectrumLoader._to_list(user_input)
        if len(input_list) != expected_length:
            raise ValueError(
                f"Length mismatch for parameter '{param_name}'. "
                f"Expected {expected_length} elements (to match 'pi'), but got {len(input_list)}."
            )
        return input_list

    @staticmethod
    def _resolve_and_validate(file_param, current_path, file_description):
        """
        Resolve the complete file path and confirm its physical existence on disk.

        Combines a file name with a common directory path string while safely preserving 
        the user's layout style (supporting relative and absolute paths). Verifies that 
        the target points to a real file on the filesystem.

        :param file_param: The relative filename or absolute path string of the target file.
        :type file_param: str
        :param current_path: The optional root directory path associated with this file index.
        :type current_path: str or None
        :param file_description: A clear label descriptor of the file type (used for error reporting).
        :type file_description: str

        :return: The fully resolved, normalized filesystem path to the file.
        :rtype: str

        :raises FileNotFoundError: If the file cannot be located at the resolved destination path.
        """
        # Strict Rule: If a common path is set, the filename cannot be absolute
        if current_path and os.path.isabs(file_param):
            expected_path = f"{current_path}/{file_param}".replace("//", "/")
        else:
            expected_path = os.path.join(current_path, file_param) if current_path else file_param
        
        # Verify it exists exactly there
        if not os.path.isfile(expected_path):
            raise FileNotFoundError(f"{file_description} not found at expected location: '{expected_path}'")
            
        return expected_path
 
    def _check_alldata_status(self, xspec_loaded_spectra):
        '''
        Check live PyXspec AllData database states to ensure it matches with self after running SpectrumLoader.xspec_load()

        Performs a check on `xspec.Spectrum` data entries. 
        Iterates through background files, RMFs, and ARFs to confirm that XSpec's status matches 
        the loader's definitions. This method guards against XSpec silently ignoring user files, 
        loading default pipeline fallbacks, or retaining cross-contamination from old sessions.

        :param xspec_loaded_spectra: A single active XSpec spectrum or a list of active spectra 
                                     currently registered inside `xspec.AllData`.
        :type xspec_loaded_spectra: xspec.Spectrum or list[xspec.Spectrum]

        :raises ValueError: If a mismatch is discovered (e.g., a file is loaded when it shouldn't 
                            be, a file is missing entirely, or file basenames do not match).
        '''
        # Get spectra back into a list if needed
        if isinstance(xspec_loaded_spectra, list):
            loaded_spectra = xspec_loaded_spectra
        else:
            loaded_spectra = [xspec_loaded_spectra]

        for i, s in enumerate(loaded_spectra):
            
            # 1. Verify the spectrum itself
            # If the spectrum file was not a proper spectrum, xspec_load will have already 
            # thrown an error
            if s.fileName != self.pi[i]:
                raise ValueError(
                    f"[wrap_xspec] Error at index {i}: Loader expected spectrum '{self.pi[i]}', "
                    f"but XSpec has spectrum '{s.fileName}' loaded."
                )

            # 2. Make a dictionary of expected values
            components_to_check = [
                {
                    "name": "background",
                    "expected_val": self.back[i],
                    "get_actual": lambda spec: spec.background.fileName,
                    "missing_err_msg": "has no background"
                },
                {
                    "name": "response (RMF)",
                    "expected_val": self.rmf[i],
                    "get_actual": lambda spec: spec.response.rmf,
                    "missing_err_msg": "No response is assigned"
                },
                {
                    "name": "ARF",
                    "expected_val": self.arf[i],
                    "get_actual": lambda spec: spec.response.arf,
                    "missing_err_msg": "Response has no Arf"
                }
            ]

            # 3. Iterate through the components
            for comp in components_to_check:
                name = comp["name"]
                expected = comp["expected_val"]
                
                if expected is None:
                    # We expect XSpec to NOT have this file loaded
                    try:
                        actual = comp["get_actual"](s)
                        raise ValueError(
                            f"[wrap_xspec] Error at index {i}: Loader expected NO {name}, "
                            f"but XSpec has an active {name} loaded: '{actual}'."
                        )
                    except Exception as e:
                        if comp["missing_err_msg"] not in str(e):
                            raise e  # Re-raise unexpected structural failures
                else:
                    # We expect XSpec to HAVE this file loaded, and it must match
                    try:
                        actual = comp["get_actual"](s)

                    except Exception as e:
                        raise ValueError(
                            f"[wrap_xspec] Error at index {i}: Loader expected {name} '{expected}', "
                            f"but XSpec reports no {name} is loaded. Check that '{expected}' is really a {name} file. "
                            #f" Underlying issue: {e}"
                        )        
                    if actual != expected:
                        raise ValueError(
                            f"[wrap_xspec] Error at index {i}: Loader expected {name} '{expected}', "
                            f"but XSpec has {name} '{actual}' loaded."
                        )

