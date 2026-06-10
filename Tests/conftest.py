import os
import pytest

# Calculate paths relative to this conftest.py file
# ==========================================
# DYNAMIC PATH CALCULATION
# ==========================================
# os.path.abspath(__file__) gets the absolute path of test_file.py
# os.path.dirname(...) gets the 'Tests' directory containing it
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

# Concat to get the absolute path to 'TestInputData' regardless of where pytest is run
TEST_DATA_ROOT = os.path.join(TESTS_DIR, "TestInputData")

# Subdirectory Targets
SAME_PATH_DIR = os.path.join(TEST_DATA_ROOT, "Test_load_from_same_path")
DIFF_PATH_DIR = os.path.join(TEST_DATA_ROOT, "Test_load_from_diff_path")

FOLDER_1 = os.path.join(DIFF_PATH_DIR, "Spectrum1")
FOLDER_2 = os.path.join(DIFF_PATH_DIR, "Spectrum2")
FOLDER_3 = os.path.join(DIFF_PATH_DIR, "Spectrum3")
FOLDER_4 = os.path.join(DIFF_PATH_DIR, "Spectrum4")


# Core test file names
PI_1, PI_2, PI_3, PI_4 = "ACISspectrum1.pi", "ACISspectrum2.pi", "ACISspectrum3.pi", "ACISspectrum4.pi"
BACK_1, BACK_2, BACK_3, BACK_4 = "ACISspectrum1_bkg.pi", "ACISspectrum2_bkg.pi", "ACISspectrum3_bkg.pi", "ACISspectrum4_bkg.pi"
RMF_1, RMF_2, RMF_3, RMF_4 = "ACISspectrum1.rmf", "ACISspectrum2.rmf", "ACISspectrum3.rmf", "ACISspectrum4.rmf"
ARF_1, ARF_2, ARF_3, ARF_4 = "ACISspectrum1.arf", "ACISspectrum2.arf", "ACISspectrum3.arf", "ACISspectrum4.arf"


# ----------------------------------------------------------------------
# 1. Configuration: Single Spectrum
# ----------------------------------------------------------------------
@pytest.fixture(scope="session")
def load_spectrum_config():
    """Configuration for a single standalone spectrum inside the same folder."""    
    return {
        "same_path_string": SAME_PATH_DIR,
        "diff_path_list": [FOLDER_1, FOLDER_2, FOLDER_3, FOLDER_4],
        "single_pi_string": PI_1,
        "single_back_string": BACK_1,
        "single_rmf_string": RMF_1,
        "single_arf_string": ARF_1,
        "pi_list": [PI_1, PI_2, PI_3, PI_4],
        "back_list": [BACK_1, BACK_2, BACK_3, BACK_4],
        "rmf_list": [RMF_1, RMF_2, RMF_3, RMF_4],
        "arf_list": [ARF_1, ARF_2, ARF_3, ARF_4],
        }



@pytest.fixture(scope="session", autouse=True)
def test_spectrum_config(load_spectrum_config):
    """Verifies that the test directory and mandatory files exist once per test session."""

    conf = load_spectrum_config

    # Check that the relevant TestInputData exists
    if not os.path.exists(conf["same_path_string"]):
        pytest.fail(f"Test directory '{conf['same_path_string']}' does not exist.")

    for i in conf["diff_path_list"]:
        if not os.path.exists(i):
            pytest.fail(f"Test directory '{i}' does not exist.")

    # Check that all of the relevant files exists

    single_keys = [key for key in conf if "single" in key]    
    for key in single_keys:    
    
        file = os.path.join(conf["same_path_string"], conf[key])
        if not os.path.isfile(file):
            pytest.fail(f"Mandatory test file '{file}' missing.")
    
    for key in ["pi_list", "back_list", "rmf_list", "arf_list"]:
        for j in range(0,len(conf[key])):  
            file = os.path.join(conf["same_path_string"], conf[key][j])
            if not os.path.isfile(file):
                pytest.fail(f"Mandatory test file '{file}' missing.")            
            file = os.path.join(conf["diff_path_list"][j], conf[key][j])
            if not os.path.isfile(file):
                pytest.fail(f"Mandatory test file '{file}' missing.")            
