transformer = {
    'inputs': ['pdb_aas', 'q8s', 'msas'],
    'targets': ['phis', 'psis'],  # 'dcalphas', 'phis', 'psis'

    'dataset_reader': {
        'root': '/home/Wayne/protein/data',
        'train_file_prefix': 'training_100_',
        'total_folds': 104,
        'params': {
            'position_emb_dim': 64,  # None if not using position embedding
            'window_size': None
        },
        'train_folds': [1, 2, 3, 4, 5, 6, 7, 8, 9],
        'valid_folds': [10],
        'batch_size': 8,
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
            'dropout': 0.
        },

        'encoder': {
            'type': 'transformer',
            'num_hidden_layers': 1,
            'num_heads': 4,
            'filter_size': 256,
            'attention_dropout': 0.,
            'relu_dropout': 0.,
            'allow_ffn_pad': True,
            'postprocess_dropout': 0.
        },

        'feedforward': {
            'hid_dim': 64,
            'num_layers': 2,
            'activation': 'relu',
            'dropout': 0.,
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
                'name': 'phis'
            },
            'psis': {
                'type': 'elementwise_linear',
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
