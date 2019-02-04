# -*- coding: utf-8 -*-
import urllib.request
import shutil
import os
import numpy as np
import csv
import gzip
from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, Activation
import logging
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from keras import backend as K

"""Main module."""
# Data parameters

SCHEME = "http://enterobase.warwick.ac.uk/schemes/Escherichia.wgMLST/profiles.list.gz"
GENE_LABEL = {"locus_tag": "ECs1205", "description": "Shiga toxin 2 subunit A"}
DATA_DIR = os.path.expanduser("~")
PROFILES = os.path.join(DATA_DIR, "profiles.gz")
TOTAL_COUNT = 20000

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

logger.info("Setting Parameters:")
# Model parameters
batch_size = 1000
nb_epoch = 12
# number of convolutional filters to use
nb_filters = 32
# size of pooling area for max pooling
pool_size = (2, 2)
# convolution kernel size
kernel_size = 3
strides = 3

def f1(y_true, y_pred):
    def recall(y_true, y_pred):
        """Recall metric.

        Only computes a batch-wise average of recall.

        Computes the recall, a metric for multi-label classification of
        how many relevant items are selected.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
        recall = true_positives / (possible_positives + K.epsilon())
        return recall

    def precision(y_true, y_pred):
        """Precision metric.

        Only computes a batch-wise average of precision.

        Computes the precision, a metric for multi-label classification of
        how many selected items are relevant.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
        precision = true_positives / (predicted_positives + K.epsilon())
        return precision
    precision = precision(y_true, y_pred)
    recall = recall(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))


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
# TODO: Getting this into classes would be good
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

# TODO: Change to use data generators - iterators; so we dont run out of memory.
labelencoder_y_1 = LabelEncoder()
all_gene_labels_enc = labelencoder_y_1.fit_transform(np.array(all_gene_labels))
nb_classes = max(all_gene_labels_enc) + 1
value_train, value_test, label_train, label_test = train_test_split(
    np.array(all_gene_values), all_gene_labels_enc, test_size=0.3, random_state=0
)
value_test = StandardScaler(copy=False).fit_transform(value_test)
value_train = StandardScaler(copy=False).fit_transform(value_train)

# Build the model
logger.info("Building the ANN")
model = Sequential()
model.add(Dense(128, activation="relu", input_dim=input_length))

#model.add(Dense(128, activation="relu"))
#model.add(Dropout(0.5))
model.add(Dense(nb_classes))
model.add(Activation("softmax"))

# TODO: Where is the info about what features are guiding the prediction.
model.compile(
    loss="sparse_categorical_crossentropy", optimizer="adadelta", metrics=[f1]
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

# TODO: Export hidden features
layer_name = "activation_1"
intermediate_layer_model = Model(
    inputs=model.input, outputs=model.get_layer(layer_name).output
)
intermediate_output = intermediate_layer_model.predict(value_test)

plt.plot(intermediate_output)
plt.show()

# TODO: These graphs should write to file and should be tidy. We should show a plot for each layer.
# Plot training & validation accuracy values
plt.plot(history.history["f1"])
plt.plot(history.history["val_f1"])
plt.title("Model accuracy")
plt.ylabel("F1 Score")
plt.xlabel("Epoch")
plt.legend(["Train", "Test"], loc="upper left")
plt.show()
#
# # Plot training & validation loss values
# plt.plot(history.history["loss"])
# plt.plot(history.history["val_loss"])
# plt.title("Model loss")
# plt.ylabel("Loss")
# plt.xlabel("Epoch")
# plt.legend(["Train", "Test"], loc="upper left")
# plt.show()
