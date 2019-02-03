# -*- coding: utf-8 -*-
import urllib.request
import shutil
import os
import numpy as np
import csv
import gzip
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
import logging
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

"""Main module."""
# Data parameters

SCHEME = "http://enterobase.warwick.ac.uk/schemes/Escherichia.wgMLST/profiles.list.gz"
GENE_LABEL = {"locus_tag": "ECs1205", "description": "Shiga toxin 2 subunit A"}
DATA_DIR = os.path.expanduser("~")
PROFILES = os.path.join(DATA_DIR, "profiles.gz")
TOTAL_COUNT = 5000

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

logger.info("Setting Parameters:")
# Model parameters
batch_size = 128
nb_epoch = 12
# number of convolutional filters to use
nb_filters = 32
# size of pooling area for max pooling
pool_size = (2, 2)
# convolution kernel size
kernel_size = 3
strides = 3


if not os.path.exists(PROFILES):
    logger.info("Profile does not exist, downloading...")
    # Fetch full profile list from EnteroBase if files doesn't exist
    with urllib.request.urlopen(SCHEME) as response, open(PROFILES, "wb") as out_file:
        shutil.copyfileobj(response, out_file)

# Load data
count = 0
all_gene_values = []
all_gene_labels = []
np.random.seed(1234)
input_length = 0

logger.info("Building data columns up to {} records...".format(TOTAL_COUNT))
for x in csv.DictReader(gzip.open(PROFILES, "rt"), dialect=csv.excel_tab):
    x.pop("ST", None)  # We don't want ST used
    gene_label = int(x.pop(GENE_LABEL["locus_tag"]))
    gene_values = np.array([int(x) for x in x.values()])
    input_length = len(x.values())
    all_gene_labels.append(gene_label)
    all_gene_values.append(gene_values)
    count += 1
    if count > TOTAL_COUNT:
        break
    elif count % 1000 == 0:
        logging.info("Loaded {} records".format(count))

logger.info("Normalising and fitting data")
labelencoder_y_1 = LabelEncoder()
all_gene_labels_enc = labelencoder_y_1.fit_transform(np.array(all_gene_labels))
nb_classes = max(all_gene_labels_enc) + 1
value_train, value_test, label_train, label_test = train_test_split(
    np.array(all_gene_values), all_gene_labels_enc, test_size=0.3, random_state=0
)
sc = StandardScaler()
value_test = sc.fit_transform(value_test)
value_train = sc.fit_transform(value_train)

# Build the model
logger.info("Building the CNN")
model = Sequential()
model.add(Dense(32, activation="relu", input_dim=input_length))

model.add(Dense(128, activation="relu"))
model.add(Dropout(0.5))
model.add(Dense(nb_classes))
model.add(Activation("softmax"))

model.compile(
    loss="sparse_categorical_crossentropy", optimizer="adadelta", metrics=["accuracy"]
)
logger.info("Summary follows")
model.summary()

logger.info("Training...")
# Train & Test
history = model.fit(
    value_train,
    label_train,
    batch_size=batch_size,
    epochs=nb_epoch,
    verbose=1,
    validation_data=(value_test, label_test),
)
logger.info("Unit complete...")

# Export hidden features
