Model 4: BiGRUs and BiLSTMs

| Author     | Affiliation         | Email               |
| ---------- | ------------------- | ------------------- |
| Arjun Srivatsa | Columbia University | a.srivatsa@columbia.edu |
| Daniel Jeong | Columbia University | daniel.jeong@columbia.edu |

# Environment
| TensorFlow Ver. | CUDA Ver. | CuDNN Ver. |
| --------------- | --------- | ---------- |
| 1.13.1          | 10.0      | 7.5        |

# Training
The scripts needed for training our proposed RNN encoder-decoder networks for predicting 3D coordinates, distance matrices (dcalpha), and torsion angles are `train_rnn_coords.py`, `train_rnn_coords_dcalpha.py`, `train_rnn_torsion.py`, and `train_rnn_torsion_dcalpha.py`. More specifically,

| Training Script | Description |
| --------------- | ----------- |
| train_rnn_coords.py | 3D Coordinate Prediction |
| train_rnn_coords_dcalpha.py | 3D-Coordinate-to-Distance-Matrix Prediction |
| train_rnn_torsion_dcalpha.py | Torsion Angle and Torsion-Angle-to-Distance-Matrix Prediction |
| train_rnn_torsion.py | Torsion Angle Prediction |

The scripts assume a particular directory structure and uses the configurations (dataset path, model output path, etc.) that are defined in `core.py`. The default directory structure is as the following:
![default_directory_structure](/model_4/description/directory_structure.png)

Once everything is in place, run the following commands:
```
python3 train_rnn_coords.py
python3 train_rnn_coords_dcalpha.py
python3 train_rnn_torsion.py
python3 train_rnn_torsion_dcalpha.py
```

# Loading and Evaluating the Trained Models
The `load_eval_rnn.py` script loads the trained coords, dcalpha, and/or torsion angle models, computes the RMSD using the validation set, and generates predictions using the test set. As with the scripts used for training, all of the settings for the script are defined in `core.py`. Depending on the specific model that you want to evaluate, the following command line options are available: `--coords`, `--coords_dcalpha`, `--torsion`, `--torsion_dcalpha`. Multiple options can be given at the same time, and in the case where distance matrix and torsion angle predictions are both available (coordinate predictions included if available), a pickle file named `predictions.pkl` containing a mapping from each test protein's CASP ID to the predicted outputs are also saved. All of the predictions (along with the distance matrix colormaps) will be saved in the specified output paths.
```
python3 load_eval_rnn.py --coords
python3 load_eval_rnn.py --coords_dcalpha
python3 load_eval_rnn.py --torsion
python3 load_eval_rnn.py --torsion_dcalpha
python3 load_eval_rnn.py --coords --coords_dcalpha --torsion_dcalpha
```
