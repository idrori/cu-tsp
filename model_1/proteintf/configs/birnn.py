birnn = {
    'inputs': ['pdb_aas', 'q8s', 'msas'],
    'targets': ['dcalphas'],  # 'dcalphas', 'phis', 'psis'

    'dataset_reader': {
        'root': '/home/Wayne/protein/data',
        'train_file_prefix': 'training_100_',
        'total_folds': 104,
        'params': {
            'position_emb_dim': 64,  # None if not using position embedding
            'window_size': None
        },
        'train_folds': list(range(1, 95)),
        'valid_folds': list(range(95, 105)),
        'batch_size': 32,
        'max_length': 384,
        'normalize_target': False,
        'shuffle': True
    },

    'model': {
        'type': 'single_seq_model',
        'window_size': None,
        'embedder': {
            'embedders': {
                'position': {
                    'type': 'pass_through',
                    'input_dim': 64
                },
                'pdb_aas': {
                    'type': 'embedding',
                    'vocab_size': 21,
                    'emb_dim': 64,
                    'project_dim': None,
                    'mask_zero': True,
                },
                'q8s': {
                    'type': 'embedding',
                    'vocab_size': 8,
                    'emb_dim': 64,
                    'project_dim': None,
                    'mask_zero': True,
                },
                'msas': {
                    'type': 'dense',
                    'input_dim': 21,
                    'output_dim': 64,
                    'batch_norm': True
                }
            },
            'merge_mode': 'sum',  # 'concat', 'sum', or 'ave'
            'project_dim': None,
            'dropout': 0.2
        },

        'encoder': {
            'type': 'rnn',
            'hid_dim': 32,
            'num_layers': 1,
            'bidirectional': True,
            'batch_norm': True,
            'residual': True,
            'dropout': 0.,
            'RNN': 'CuDNNGRU'
        },

        'feedforward': {
            'hid_dim': 64,
            'num_layers': 2,
            'activation': 'relu',
            'dropout': 0.2,
            'batch_norm': True,
            'residual': True
        },

        'decoder': {
            'dcalphas': {
                'type': 'pairwise_bilinear',
                'name': 'dcalphas'
            },
            'phis': {
                'type': 'elementwise_linear',
                # 'feedforward': {
                #     'hid_dim': 64,
                #     'num_layers': 2,
                #     'activation': 'relu',
                #     'dropout': 0.2,
                #     'batch_norm': True,
                #     'residual': True
                # },
                'name': 'phis'
            },
            'psis': {
                'type': 'elementwise_linear',
                # 'feedforward': {
                #     'hid_dim': 64,
                #     'num_layers': 2,
                #     'activation': 'relu',
                #     'dropout': 0.2,
                #     'batch_norm': True,
                #     'residual': True
                # },
                'name': 'psis'
            }
        },

        'loss': {
            'dcalphas': {
                'type': 'rmse',
                'dim': 3
            },
            'phis': {
                'type': 'rmse',
                'dim': 2
            },
            'psis': {
                'type': 'rmse',
                'dim': 2
            }
        },

        # 'metrics': {
        #     'dcalphas': {
        #         'type': 'rmse',
        #         'dim': 3
        #     },
        #     'phis': {
        #         'type': 'rmse',
        #         'dim': 2
        #     },
        #     'psis': {
        #         'type': 'rmse',
        #         'dim': 2
        #     }
        # },

        'optimizer': {
            'type': 'Adam',
            'lr': 0.001
        }
    },

    'num_epochs': 1000,
    'callbacks': {
        'model_checkpoint': {
            'monitor': 'val_loss',
            'mode': 'min',
            'save_best_only': True,
            'verbose': 1
        },
        'early_stopping': {
            'monitor': 'val_loss',
            'mode': 'min',
            'patience': 20,
            'verbose': 1
        },
        'tensorboard': {
            # 'histogram_freq': 1,
            # 'write_grads': True,
            # 'write_images': True
        }
        # 'noam_lr_scheduler': {
        #     'model_size': 64,
        #     'warmup_steps': 8000,
        #     'factor': 1.,
        #     'verbose': 1
        # }
    }
}
