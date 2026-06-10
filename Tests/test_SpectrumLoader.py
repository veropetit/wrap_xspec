import xspec
xspec.Xset.allowPrompting = False
from wrap_xspec import SpectrumLoader
import pytest
import os
import numpy as np

class Test_SpectrumLoader_validate_input_type:

    def test_accept_valid_type(self):
        """Verify that the type validation accepts lists and strings"""
        input = "string"
        SpectrumLoader._validate_input_type(input, 'param', allow_none=False)
        SpectrumLoader._validate_input_type(input, 'param', allow_none=True)
        input = None
        SpectrumLoader._validate_input_type(input, 'param', allow_none=True)
        input = ("string1", "string2")
        SpectrumLoader._validate_input_type(input, 'param', allow_none=False)
        SpectrumLoader._validate_input_type(input, 'param', allow_none=True)
        input = ("string1", None)
        SpectrumLoader._validate_input_type(input, 'param', allow_none=True)
        input = ["string1", "string2"]
        SpectrumLoader._validate_input_type(input, 'param', allow_none=False)
        SpectrumLoader._validate_input_type(input, 'param', allow_none=True)
        input = ["string1", None]
        SpectrumLoader._validate_input_type(input, 'param', allow_none=True)
        input = np.array(["string1", "string2"])
        SpectrumLoader._validate_input_type(input, 'param', allow_none=False)
        SpectrumLoader._validate_input_type(input, 'param', allow_none=True)
        input = np.array(["string1", None])
        SpectrumLoader._validate_input_type(input, 'param', allow_none=True)

    def test_rejects_invalid_base_types(self):
        """Verify that the type validation crashes if invalid wrapper structures are passed."""
        # Dictionaries should fail
        with pytest.raises(TypeError, match="must be a string or a sequence of strings"):
            input = pi={"file": "spec.pi"}
            SpectrumLoader._validate_input_type(input, 'param', allow_none=True)

        # Integers should fail
        with pytest.raises(TypeError, match="must be a string or a sequence of strings"):
            input = 123
            SpectrumLoader._validate_input_type(input, 'param', allow_none=True)

        # A None when not allowed
        with pytest.raises(TypeError, match="Parameter 'param' cannot be None."):
            SpectrumLoader._validate_input_type(None, 'param', allow_none=False)

    def test_rejects_invalid_sequence_elements(self):
        """Verify a sequence fails if it contains nested non-string objects."""
        # List containing an integer
        with pytest.raises(TypeError, match="All elements in parameter 'param' must be strings. Found type 'int' instead."):
            SpectrumLoader._validate_input_type(["spec1.pi", 42, "spec3.pi"], 'param', allow_none=True)

        # Numpy array containing floats
        with pytest.raises(TypeError, match="All elements in parameter 'param' must be strings. Found type 'float64' instead."):
            SpectrumLoader._validate_input_type(np.array([1.0, 2.5, 3.1]), 'param', allow_none=True)

        # A list containing a None when not supposed to. 
        with pytest.raises(TypeError, match="All elements in parameter 'param' must be strings. Found type 'NoneType' instead."):
            SpectrumLoader._validate_input_type(["spec1.pi", None, "spec3.pi"], 'param', allow_none=False)

class Test_SpectrumLoader_to_list:
    @pytest.mark.parametrize('input, output', [
        ('a', ['a']),
        (None, [None]),
        (['a', 'b'], ['a', 'b']),
        (('a','b'), ['a', 'b']),
        (np.array(['a', 'b']), ['a','b']),
    ])
    def test_to_list_sucess(self, input, output):
        '''Def the to_list static method'''
        # The input is already valiated upstream by _validate_input_type
        result = SpectrumLoader._to_list(input)
        assert result == output
        
class Test_broadcast_input:
    @pytest.mark.parametrize('input, output, exp_len', [
        ('a', ['a', 'a', 'a'], 3),
        (None, [None, None, None], 3),
        (['a', 'b'], ['a', 'b'], 2),
        (('a','b'), ['a', 'b'], 2),
        (np.array(['a', 'b']), ['a','b'], 2),
    ])
    def test_broadcast_input_success(self, input, output, exp_len):
        result = SpectrumLoader._broadcast_input(input, exp_len, "param")
        assert result == output     

    def test_broadcast_input_failure(self):
        with pytest.raises(ValueError, match=r"Length mismatch for parameter 'param'. Expected 5 elements \(to match 'pi'\), but got 3."):
            SpectrumLoader._broadcast_input(['a', 'b', 'c'], 5, "param")

class Test_resolve_and_validate():

    def test_preserves_relative_path_with_common_path(self, load_spectrum_config):
        """If a relative common path and filename are given, the combined relative path is retained."""
        conf = load_spectrum_config
        
        # Force a relative path styling
        relative_path = os.path.relpath(conf["same_path_string"])
        
        result = SpectrumLoader._resolve_and_validate(
            conf["single_pi_string"], relative_path, "param"
        )
        
        expected_combined = os.path.join(relative_path, conf["single_pi_string"])
        assert result == expected_combined # loader forces lists
        # Ensure it didn't accidentally convert to an absolute path string
        assert not result.startswith("/")

    def test_preserves_absolute_path_without_common_path(self, load_spectrum_config):
        """If no common path is given and the file is absolute, the absolute string is retained."""
        conf = load_spectrum_config

        absolute_pi_path = os.path.abspath(os.path.join(conf["same_path_string"], conf["single_pi_string"]))
        
        result = SpectrumLoader._resolve_and_validate(
            absolute_pi_path, None, "param"
        )
                        
        assert result == absolute_pi_path
        assert result.startswith("/")   

    def test_strict_failure_on_mixed_locations(self, load_spectrum_config):
        """
        If a user passes an absolute file path but ALSO specifies a common path,
        the code should strictly concatenate them, fail to find the file, and raise an error.
        """
        conf = load_spectrum_config

        absolute_pi_path = os.path.abspath(os.path.join(conf["same_path_string"], conf["single_pi_string"]))
        fake_common_path = conf["same_path_string"]
        
        # This should fail because it expects the file at 'SomeCommonDirectory//absolute/path/to/file.pi'
        with pytest.raises(FileNotFoundError) as exc_info:
            SpectrumLoader(path=fake_common_path, pi=absolute_pi_path)
            
        # Verify the error message reflects the strict, non-second-guessed path
        expected_broken_path = os.path.join(fake_common_path, absolute_pi_path)
        assert expected_broken_path in str(exc_info.value)

class Test_SpectrumLoader_Init:

    def test_integration_accepts_valid_types_and_broadcast(self, monkeypatch):
        """Verify that strings, lists, tuples, and numpy arrays of strings pass type validation."""
        # Mock _resolve_and_validate to always pass so we can focus solely on type checking
        monkeypatch.setattr(SpectrumLoader, "_resolve_and_validate", lambda self, f, p, d: f)
        monkeypatch.setattr(os.path, "exists", lambda path: True)

        # 1. Test standard string
        loader = SpectrumLoader(pi="spec.pi")
        assert loader._pi_list == ["spec.pi"]

        # 2. Test Tuple input
        loader = SpectrumLoader(pi=("spec1.pi", "spec2.pi"))
        assert loader._pi_list == ["spec1.pi", "spec2.pi"]

        # 3. Test NumPy array input
        np_input = np.array(["spec1.pi", "spec2.pi", "spec3.pi"])
        loader = SpectrumLoader(pi=np_input)
        assert loader._pi_list == ["spec1.pi", "spec2.pi", "spec3.pi"]

    def test_integration_failed_types(self):
        """Verify initialization immediately crashes if invalid wrapper structures are passed.
        and test that the allow_none are set correctly"""
        
        with pytest.raises(TypeError, match="must be a string or a sequence of strings"):
            SpectrumLoader(pi=12345)

    def test_successful_init_1spec_mandatory_only(self, load_spectrum_config):
        """Test initialization with only the mandatory path and PI file."""
        conf = load_spectrum_config
        loader = SpectrumLoader(path=conf["same_path_string"], pi=conf["single_pi_string"])
        
        assert loader._path_list == [conf["same_path_string"]]
        assert loader.pi == [os.path.join(conf["same_path_string"], conf["single_pi_string"])]
        assert loader.back == [None]
        assert loader.rmf == [None]
        assert loader.arf == [None]

    def test_successful_init_1spec_all_files(self, load_spectrum_config):
        """Test initialization when all optional files are provided and valid."""
        conf = load_spectrum_config

        loader = SpectrumLoader(
            path=conf["same_path_string"],
            pi=conf["single_pi_string"],
            back=conf["single_back_string"],
            rmf=conf["single_rmf_string"],
            arf=conf["single_arf_string"]
        )
        
        assert loader._path_list == [conf["same_path_string"]]
        assert loader.pi == [os.path.join(conf["same_path_string"], conf["single_pi_string"])]
        assert loader.back == [os.path.join(conf["same_path_string"], conf["single_back_string"])]
        assert loader.rmf == [os.path.join(conf["same_path_string"], conf["single_rmf_string"])]
        assert loader.arf == [os.path.join(conf["same_path_string"], conf["single_arf_string"])]

    def test_successful_init_4spec_all_files_one_path(self, load_spectrum_config):
        """Test initialization when all optional files are provided and valid."""
        conf = load_spectrum_config

        loader = SpectrumLoader(
            path=conf["same_path_string"],
            pi=conf["pi_list"],
            back=conf["back_list"],
            rmf=conf["rmf_list"],
            arf=conf["arf_list"]
        )
        
        assert loader._path_list == [conf["same_path_string"]]*len(conf["pi_list"])
        assert loader.pi == [ os.path.join(conf["same_path_string"],x) for x in conf["pi_list"]]
        assert loader.back == [os.path.join(conf["same_path_string"], x) for x in conf["back_list"]]
        assert loader.rmf == [os.path.join(conf["same_path_string"],x) for x in conf["rmf_list"]]
        assert loader.arf == [os.path.join(conf["same_path_string"],x) for x in conf["arf_list"]]

    def test_successful_init_4spec_mandatory_only_one_path(self, load_spectrum_config):
        """Test initialization when all optional files are provided and valid."""
        conf = load_spectrum_config

        loader = SpectrumLoader(
            path=conf["same_path_string"],
            pi=conf["pi_list"],
            rmf=conf["single_rmf_string"]
        )
        
        assert loader._path_list == [conf["same_path_string"]]*len(conf["pi_list"])
        assert loader.pi == [ os.path.join(conf["same_path_string"],x) for x in conf["pi_list"]]
        assert loader.back == [None]*len(conf["pi_list"])
        assert loader.rmf == [os.path.join(conf["same_path_string"],conf['single_rmf_string'])]*len(conf["pi_list"])
        assert loader.arf == [None]*len(conf["pi_list"])

    def test_successful_init_4spec_all_files_4_path(self, load_spectrum_config):
        """Test initialization when all optional files are provided and valid."""
        conf = load_spectrum_config

        loader = SpectrumLoader(
            path=conf["diff_path_list"],
            pi=conf["pi_list"],
            back=conf["back_list"],
            rmf=conf["rmf_list"],
            arf=conf["arf_list"]
        )
        
        assert loader._path_list == conf["diff_path_list"]
        assert loader.pi == [ os.path.join(y,x) for y,x in zip(conf['diff_path_list'],conf["pi_list"])]
        assert loader.back == [ os.path.join(y,x) for y,x in zip(conf['diff_path_list'],conf["back_list"])]
        assert loader.rmf == [ os.path.join(y,x) for y,x in zip(conf['diff_path_list'],conf["rmf_list"])]
        assert loader.arf == [ os.path.join(y,x) for y,x in zip(conf['diff_path_list'],conf["arf_list"])]

    def test_type_validate_failure_integration(self, load_spectrum_config):
        conf = load_spectrum_config

        # the static method _type_validate has been tested on its own alread
        # here we are just testing the integration
        # the allowed None were already tested also in the 'mandatory_only above
        with pytest.raises(TypeError):
            SpectrumLoader(pi=None)
        with pytest.raises(TypeError):
            SpectrumLoader(pi=123)
        with pytest.raises(TypeError):
            SpectrumLoader(pi=conf['single_pi_string'], path=conf['same_path_string'], arf=123)
        with pytest.raises(TypeError):
            SpectrumLoader(pi=conf['single_pi_string'], path=conf['same_path_string'], rmf=123)
        with pytest.raises(TypeError):
            SpectrumLoader(pi=conf['single_pi_string'], path=conf['same_path_string'], back=123)
        with pytest.raises(TypeError):
            SpectrumLoader(pi=conf['single_pi_string'], path=123)

    def test_integration_invalid_common_path_raises_error(self):
        """Test that a non-existent directory path raises FileNotFoundError."""
        # the validity of a single common directory is checked first
        fake_path = "this_directory_does_not_exist_123"
        with pytest.raises(FileNotFoundError) as exc_info:
            SpectrumLoader(path=fake_path, pi='FakeFile')
        assert "path directory does not exist" in str(exc_info.value)

    def test_missing_pi_file_raises_error(self, load_spectrum_config):
        """Test that a missing mandatory PI file raises FileNotFoundError."""
        # This test doesn't need setup_test_environment because it expects a failure
        fake_pi = "non_existent_spectrum.pi"
        
        with pytest.raises(FileNotFoundError) as exc_info:
            SpectrumLoader(path=load_spectrum_config["same_path_string"], pi=fake_pi)
            
        expected_target = os.path.join(load_spectrum_config["same_path_string"], fake_pi)
        assert f"Mandatory spectrum file not found at expected location: '{expected_target}'" in str(exc_info.value)

    @pytest.mark.parametrize("missing_file_kwargs, error_message_snippet", [
        ({"back": "missing_back.bak"}, "Background file not found"),
        ({"rmf": "missing_response.rmf"}, "RMF file not found"),
        ({"arf": "missing_ancillary.arf"}, "ARF file not found")
    ])
    def test_missing_optional_files_raise_error(self, load_spectrum_config, missing_file_kwargs, error_message_snippet):
        """
        Test each optional file individually to ensure it raises FileNotFoundError 
        if specified but missing.
        """
        conf = load_spectrum_config
        # Merge the mandatory argument with the parameterized missing optional file
        kwargs = {"path": conf['same_path_string'], "pi": conf['single_pi_string']}
        kwargs.update(missing_file_kwargs)
        
        with pytest.raises(FileNotFoundError) as exc_info:
            SpectrumLoader(**kwargs)
            
        assert error_message_snippet in str(exc_info.value)

class Test_xspec_load:

    def test_pass(self):
        pass

    # def test_load_success_with_files_in_directory(self, monkeypatch, load_spectrum_config):
    #     '''Test opening with mandatory keywords, 
    #     and that the back, arf, and responses automatically loaded'''
    #     xspec.AllData.clear()
    #     # Moving to path
    #     starting_path = os.getcwd()

    #     monkeypatch.chdir(load_spectrum_config['path'])

    #     loader = SpectrumLoader(
    #         pi=load_spectrum_config["pi"],
    #     )
    
    #     s = loader.xspec_load()   

    #     assert s.fileName == load_spectrum_config["pi"] 

    #     # These have been automatically loaded. 
    #     #assert s.response.rmf == load_spectrum_config['rmf']
    #     #assert s.response.arf == load_spectrum_config['arf']
    #     #assert s.background.fileName == load_spectrum_config['back']

    # def test_load_success_with_all_files(self, load_spectrum_config):
    #     '''Test opening with all keywords, '''
    #     xspec.AllData.clear()
    #     starting_path = os.getcwd()
    #     loader = SpectrumLoader(
    #         path=load_spectrum_config["path"],
    #         pi=load_spectrum_config["pi"],
    #         rmf=load_spectrum_config["rmf"],
    #         arf=load_spectrum_config["arf"],
    #         back=load_spectrum_config["back"]
    #     )
    #     s = loader.xspec_load()

    #     # Check that we are back in the correct directory. 
    #     assert starting_path == os.getcwd()    

    #     assert s.fileName == load_spectrum_config["pi"] 

    #     # These have been automatically loaded. 
    #     assert s.response.rmf == load_spectrum_config['rmf']
    #     assert s.response.arf == load_spectrum_config['arf']
    #     assert s.background.fileName == load_spectrum_config['back']

    # def test_easy_groups(self, load_spectrum_config):
    #     xspec.AllData.clear()
    #     loader = SpectrumLoader(
    #         path=load_spectrum_config["path"],
    #         pi=load_spectrum_config["pi"],
    #     )
    #     s1 = loader.xspec_load()
    #     s2 = loader.xspec_load()

    #     # spectra should be in datagroup 1
    #     assert s1.dataGroup == 1
    #     assert s1.index == 1
    #     assert s2.dataGroup == 1
    #     assert s2.index == 2

    #     assert xspec.AllData.nGroups == 1
    #     assert xspec.AllData.nSpectra == 2

    # def test_load_success_with_no_back(self, load_spectrum_config):
    #     xspec.AllData.clear()
    #     loader = SpectrumLoader(
    #         path=load_spectrum_config["path"],
    #         pi=load_spectrum_config["pi"],
    #         back=''
    #     )
    #     s = loader.xspec_load()  
    #     with pytest.raises(Exception) as exc_info:
    #         # If I try to access the background, it will fail
    #         print(s.background.fileName)
    #     assert str(exc_info.value) == "Error: Spectrum has no background."

    # def test_load_success_with_no_arf(self, load_spectrum_config):
    #     xspec.AllData.clear()
    #     loader = SpectrumLoader(
    #         path=load_spectrum_config["path"],
    #         pi=load_spectrum_config["pi"],
    #         arf=''
    #     )
    #     s = loader.xspec_load()  
        
    #     with pytest.raises(Exception) as exc_info:
    #         # If I try to access the arf, it will fail
    #         print(s.response.arf)
    #     assert str(exc_info.value) == "Error: Response has no Arf."

    # def test_load_success_with_no_rmf(self, load_spectrum_config):
    #     xspec.AllData.clear()
    #     loader = SpectrumLoader(
    #         path=load_spectrum_config["path"],
    #         pi=load_spectrum_config["pi"],
    #         rmf=''
    #     )
    #     s = loader.xspec_load()  
        
    #     with pytest.raises(Exception) as exc_info:
    #         # If I try to access the arf, it will fail
    #         print(s.response.rmf)
    #     assert str(exc_info.value) == "Error: No response is assigned to source 1 for spectrum 1"

    # def test_load_success_with_wrong_header_files(self, load_spectrum_config):
    #     xspec.AllData.clear()
    #     xspec.Xset.allowPrompting = False
    #     loader = SpectrumLoader(
    #         path=load_spectrum_config["path"],
    #         pi='ChandraACIS_withWrongAnciliaryFileKeywords.pi',
    #         rmf=load_spectrum_config["rmf"],
    #         arf=load_spectrum_config["arf"],
    #         back=load_spectrum_config["back"]
    #     )
    #     s = loader.xspec_load()

class Test_Check_AllData_Status:

    def test_success(self, load_spectrum_config):
        xspec.AllData.clear()
        loader = SpectrumLoader(
            path=load_spectrum_config["same_path_string"],
            pi=load_spectrum_config["single_pi_string"],
            rmf=load_spectrum_config["single_rmf_string"],
            arf=load_spectrum_config["single_arf_string"],
            back=load_spectrum_config["single_back_string"]
        )
        s = xspec.Spectrum(dataFile=loader.pi[0], respFile=loader.rmf[0], arfFile=loader.arf[0], backFile=loader.back[0])
        
        loader._check_alldata_status(s)       

    def test_mismatched_pi(self, load_spectrum_config):
        xspec.AllData.clear()
        loader = SpectrumLoader(
            path=load_spectrum_config["same_path_string"],
            pi=load_spectrum_config["single_pi_string"],
            rmf=load_spectrum_config["single_rmf_string"],
            arf=load_spectrum_config["single_arf_string"],
            back=load_spectrum_config["single_back_string"]
        )
        s = xspec.Spectrum(dataFile=loader.pi[0], respFile=loader.rmf[0], arfFile=loader.arf[0], backFile=loader.back[0])
        loader.pi = ['AnotherPi']

        with pytest.raises(ValueError) as excinfo:
            loader._check_alldata_status(s) 
        assert "Error at index 0: Loader expected spectrum" in str(excinfo.value)
        assert "but XSpec has spectrum" in str(excinfo.value)

    @pytest.mark.parametrize(
        "loader_back, xspec_back, expected_snippets",
        [
            # Case 1: Loader has None, but XSpec has a real background
            (None, 
             "single_back_string", 
             ["Error at index 0: Loader expected NO background, but XSpec has an active background loaded"]),
            
            # Case 2: Loader has a real background, but XSpec has None
            ("single_back_string", 
             None, 
             ["Error at index 0: Loader expected background", "but XSpec reports no background is loaded"]),
            
            # Case 3: Both have background, but they do not match
            ("AnotherBack", 
             "single_back_string", 
             ["Error at index 0: Loader expected background", "but XSpec has background"]),
        ]
    )
    def test_mismatched_background(self, load_spectrum_config, loader_back, xspec_back, expected_snippets):
        '''Test of the _check_alldata_status. Because this helper function is called in SpectrumLoader.xspec_load(), 
        I cannot use xspec_load to populate xspec. So I use Spectrum directly'''
        
        xspec.AllData.clear()
        
        # SpectrumLoader with all of the correct files
        # so that it resolves all the file names/paths correctly. 
        base_args = {
            "path": load_spectrum_config["same_path_string"], 
            "pi": load_spectrum_config["single_pi_string"],
            "rmf": load_spectrum_config["single_rmf_string"], 
            "arf": load_spectrum_config["single_arf_string"],
            "back": load_spectrum_config["single_back_string"]
        }
        loader = SpectrumLoader(**base_args)

        # Use the resolved paths to load xspec directly with the xspec.Spectrum
        x_back = os.path.join(base_args['path'],load_spectrum_config.get(xspec_back, None)) if isinstance(xspec_back, str) else None
        s = xspec.Spectrum(dataFile=loader.pi[0], respFile=loader.rmf[0], arfFile=loader.arf[0], backFile=x_back)

        # Manually change the background in l_back
        # (so it won't trigger a Filenotfound)
        l_back = load_spectrum_config.get(loader_back, loader_back) # Resolves string key or fallback to direct value
        loader.back = np.atleast_1d(l_back)

        # 4. Assert the validation failure behavior
        with pytest.raises(ValueError) as excinfo:
            loader._check_alldata_status(s)
            
        error_msg = str(excinfo.value)
        for snippet in expected_snippets:
            assert snippet in error_msg

    @pytest.mark.parametrize(
        "loader_rmf, xspec_rmf, expected_snippets",
        [
            # Case 1: Loader has None, but XSpec has a real response
            (None, 
             "single_rmf_string", 
             ["Error at index 0: Loader expected NO response (RMF), but XSpec has an active response (RMF) loaded"]),
            
            # Case 2: Loader has a real response, but XSpec has None
            ("single_rmf_string", 
             None, 
             ["Error at index 0: Loader expected response (RMF)", "but XSpec reports no response (RMF) is loaded"]),
            
            # Case 3: Both have response, but they do not match
            ("AnotherRMF", 
             "single_rmf_string", 
             ["Error at index 0: Loader expected response (RMF)", "but XSpec has response (RMF)"]),
        ]
    )
    def test_mismatched_response(self, load_spectrum_config, loader_rmf, xspec_rmf, expected_snippets):
        '''Test of the _check_alldata_status. Because this helper function is called in SpectrumLoader.xspec_load(), 
        I cannot use xspec_load to populate xspec. So I use Spectrum directly'''
        
        xspec.AllData.clear()
        
        # SpectrumLoader with all of the correct files
        # so that it resolves all the file names/paths correctly. 
        base_args = {
            "path": load_spectrum_config["same_path_string"], 
            "pi": load_spectrum_config["single_pi_string"],
            "rmf": load_spectrum_config["single_rmf_string"], 
            "arf": None,
            "back": load_spectrum_config["single_back_string"]

        }
        loader = SpectrumLoader(**base_args)

        # Use the resolved paths to load xspec directly with the xspec.Spectrum
        x_rmf = os.path.join(base_args['path'],load_spectrum_config.get(xspec_rmf, None)) if isinstance(xspec_rmf, str) else None
        s = xspec.Spectrum(dataFile=loader.pi[0], respFile=x_rmf, arfFile=loader.arf[0], backFile=loader.back[0])

        # Manually change the response in l_back
        # (so it won't trigger a Filenotfound)
        l_rmf = load_spectrum_config.get(loader_rmf, loader_rmf) # Resolves string key or fallback to direct value
        loader.rmf = np.atleast_1d(l_rmf)

        # 4. Assert the validation failure behavior
        with pytest.raises(ValueError) as excinfo:
            loader._check_alldata_status(s)
            
        error_msg = str(excinfo.value)
        for snippet in expected_snippets:
            assert snippet in error_msg

    @pytest.mark.parametrize(
        "loader_arf, xspec_arf, expected_snippets",
        [
            # Case 1: Loader has None, but XSpec has a real arf
            (None, 
             "single_arf_string", 
             ["Error at index 0: Loader expected NO ARF, but XSpec has an active ARF loaded"]),
            
            # Case 2: Loader has a real response, but XSpec has None
            ("single_arf_string", 
             None, 
             ["Error at index 0: Loader expected ARF", "but XSpec reports no ARF is loaded"]),
            
            # Case 3: Both have arf, but they do not match
            ("AnotherARF", 
             "single_arf_string", 
             ["Error at index 0: Loader expected ARF", "but XSpec has ARF"]),
        ]
    )
    def test_mismatched_arf(self, load_spectrum_config, loader_arf, xspec_arf, expected_snippets):
        '''Test of the _check_alldata_status. Because this helper function is called in SpectrumLoader.xspec_load(), 
        I cannot use xspec_load to populate xspec. So I use Spectrum directly'''
        
        xspec.AllData.clear()
        
        # SpectrumLoader with all of the correct files
        # so that it resolves all the file names/paths correctly. 
        base_args = {
            "path": load_spectrum_config["same_path_string"], 
            "pi": load_spectrum_config["single_pi_string"],
            "rmf": load_spectrum_config["single_rmf_string"], 
            "arf": load_spectrum_config["single_arf_string"],
            "back": load_spectrum_config["single_back_string"]
        }
        loader = SpectrumLoader(**base_args)

        # Use the resolved paths to load xspec directly with the xspec.Spectrum
        x_arf = os.path.join(base_args['path'],load_spectrum_config.get(xspec_arf, None)) if isinstance(xspec_arf, str) else None
        s = xspec.Spectrum(dataFile=loader.pi[0], respFile=loader.rmf[0], arfFile=x_arf, backFile=loader.back[0])

        # Manually change the response in l_back
        # (so it won't trigger a Filenotfound)
        l_arf = load_spectrum_config.get(loader_arf, loader_arf) # Resolves string key or fallback to direct value
        loader.arf = np.atleast_1d(l_arf)

        # 4. Assert the validation failure behavior
        with pytest.raises(ValueError) as excinfo:
            loader._check_alldata_status(s)
            
        error_msg = str(excinfo.value)
        for snippet in expected_snippets:
            assert snippet in error_msg

    def test_integration_with_loader(self, load_spectrum_config):
        '''Test the integration of _check_alldata_status into xspec_load'''
        xspec.AllData.clear()
        
        # SpectrumLoader with all but one of the correct files 
        base_args = {
            "path": load_spectrum_config["same_path_string"], 
            "pi": load_spectrum_config["single_pi_string"],
            "rmf": load_spectrum_config["single_rmf_string"], 
            "arf": load_spectrum_config["single_arf_string"],
            "back": load_spectrum_config["single_arf_string"] # <- back is in fact a arf file
        }
        loader = SpectrumLoader(**base_args)      

        with pytest.raises(ValueError) as excinfo:
            s = loader.xspec_load()  
        assert "Error at index 0: Loader expected background" in str(excinfo.value)
        assert "but XSpec reports no background is loaded" in str(excinfo.value)