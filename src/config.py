import os
###################################
# Directory setup
###################################

HOME_PATH = os.environ["PROTONDDR"]
DATA_PATH = HOME_PATH + "/repos/MUSICian_v2/musician/data"
RAW_DATA_PATH = DATA_PATH + "/raw/{}"
INTERIM_DATA_PATH = DATA_PATH + "/interim/"
PROCESSED_DATA_PATH = DATA_PATH + "/processed/{}"
EXTERNAL_DATA_PATH = DATA_PATH + "/external/"
HUSSMANN_PATH = EXTERNAL_DATA_PATH + "Hussmann2021/Supplementary_Materials/1-s2.0-S0092867421011764-mmc{}.xlsx"
CUSTOM_DB_FILE = "MBCrisprMBAgain_1.2_cseale.db"
ORIGINAL_DB_FILE = "MBCrisprMBAgain_1.2.db"
VALIDATION_DB_FILE = "MBCrisprMBSubscreeen_Full_2.db"
TUD_VALIDATION_DB_FILE = "MBCrisprMBSubscreeen_TUD.db"
ADAMSON_DB_FILE = "Adamson.db"
RAW_FILE = "MB{}_reducedCols.txt"


###################################
# Filtering parameters
###################################

CORR_THRESHOLD = .6
BAD_TARGET_THRESHOLD = 2 
NUM_GRNAS_THRESHOLD = 2

###################################
# General parameters
###################################

# FDR_ALPHA = 0.1
# NUM_OUTLYING_SAMPLES = 1

###################################
# Experiment configuations 
# Uncomment as needed
###################################

SAMPLE = "sample"
TARGET = "target"
TARGETDIFF = "targetdiff"
PAIRED_SAMPLES = "paired-replicates"

# replicate profiles with robust cov outlier detection
# PROFILE_TYPE = SAMPLE
# FILTER_COUNT = 700
# METHOD = "pca"

# PROFILE_TYPE = SAMPLE
# FILTER_COUNT = 700
# METHOD = "robust_cov"

# PROFILE_TYPE = TARGET
# FILTER_COUNT = 700
# METHOD = "robust_cov"

PROFILE_TYPE = PAIRED_SAMPLES
FILTER_COUNT = 700
METHOD = "robust_cov"


def get_raw_file_for_sample(sample):
    return RAW_DATA_PATH.format(RAW_FILE.format(sample))

def get_db_file(use_custom=True):
    if use_custom:
        return PROCESSED_DATA_PATH.format(CUSTOM_DB_FILE)
    else:
        return PROCESSED_DATA_PATH.format(ORIGINAL_DB_FILE)

def get_validation_db_file(experiment=None):
    if experiment is None:
        return PROCESSED_DATA_PATH.format(VALIDATION_DB_FILE)
    elif experiment == "TUD":
        return PROCESSED_DATA_PATH.format(TUD_VALIDATION_DB_FILE)

def get_adamson_db_file():
    return PROCESSED_DATA_PATH.format(ADAMSON_DB_FILE)

def get_hussmann_supplementary_xlsx(number):
    return HUSSMANN_PATH.format(number)

def get_interim_dir():
    return INTERIM_DATA_PATH

def get_processed_dir():
    return PROCESSED_DATA_PATH

def get_common_barcodes():
    return HOME_PATH + "/repos/MUSICian_v2/musician/src/data/filtered_barcodes.{}.tsv".format(FILTER_COUNT)

def get_validation_barcodes(downsample=None):
    if downsample is None:
        downsample = "full"
    return HOME_PATH + "/repos/MUSICian_v2/musician/src/data/filtered_validation_barcodes.{}.{}.tsv".format(FILTER_COUNT, downsample)

def get_adamson_barcodes():
    return HOME_PATH + "/repos/MUSICian_v2/musician/src/data/filtered_adamson_barcodes.tsv"

def get_adamson():
    return HOME_PATH + "/repos/MUSICian_v2/musician/src/data/adamson.txt"

def get_experiment_name(profile_type = None, method = None, filter_count = None):
    return "{}.{}.{}".format(PROFILE_TYPE if profile_type is None else profile_type,
        METHOD if method is None else method,
        FILTER_COUNT if filter_count is None else filter_count)

def get_experiment_artifacts(profile_type = None, method = None, filter_count = None):
    folder = HOME_PATH + "/repos/MUSICian_v2/musician/notebooks/exploratory/outlier_detection/artifacts/{}/".format(get_experiment_name(profile_type, method, filter_count))
    if not os.path.exists(folder):
        os.mkdir(folder)
    return HOME_PATH + "/repos/MUSICian_v2/musician/notebooks/exploratory/outlier_detection/artifacts/{}/".format(get_experiment_name(profile_type, method, filter_count))

