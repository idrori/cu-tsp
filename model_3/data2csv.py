import pickle, sys, os
import pandas as pd
import numpy as np

# load data
model_dir = sys.argv[1]  # or your model dir 
datafile = os.path.join(sys.argv[1], 'predictions.pkl')
with open(datafile, "rb") as f:
    data = pickle.load(f)
dist, psi, phi = data

# load protein names
testfile = "./data/test.csv"
test_input = pd.read_csv(testfile, header=None)
protein_names = np.array(test_input.iloc[:,1])
protein_len = np.array(test_input.iloc[:,2])

# concatenate all output to one-dimensional
all_data = []
all_names = []
for i, pname in enumerate(protein_names):
    dist_flat = dist[i].ravel()
    array = np.concatenate([dist_flat, psi[i], phi[i]])
    all_data.append(array)

    length = protein_len[i]
    dist_names = ["{}_d_{}_{}".format(pname, i + 1, j + 1) for i in range(length) for
            j in range(length)]

    psi_names = ["{}_psi_{}".format(pname, i + 1) for i in range(length)]
    phi_names = ["{}_phi_{}".format(pname, i + 1) for i in range(length)]
    row_names = np.array(dist_names + psi_names + phi_names)
    all_names.append(row_names)

all_data = np.concatenate(all_data)
all_names = np.concatenate(all_names)
output = {"Id": all_names, "Predicted": all_data}
output = pd.DataFrame(output)
output.to_csv(os.path.join(model_dir, "submission.csv"), index=False)



