from tensorflow.keras import layers


class FeedForward:
    @classmethod
    def from_params(cls, params):
        input_dim = params['input_dim']
        hid_dim = params['hid_dim']
        num_layers = params.get('num_layers', 1)
        activation = params.get('activation', 'relu')
        dropout = params.get('dropout', 0.)
        batch_norm = params.get('batch_norm', False)
        residual = params.get('residual', False)

        if residual:
            hid_dim = input_dim
        if isinstance(hid_dim, int):
            hid_dim = [hid_dim] * num_layers
        if isinstance(activation, str):
            activation = [activation] * num_layers
        assert (isinstance(hid_dim, list) and len(hid_dim) == num_layers)
        assert (isinstance(activation, list) and len(activation) == num_layers)

        def _feedforward(x):
            in_dim = input_dim
            for out_dim, sigma in zip(hid_dim, activation):
                y = x
                if batch_norm:
                    x = layers.BatchNormalization(axis=-1)(x)
                x = layers.Dense(out_dim, input_shape=(in_dim,))(x)
                x = layers.Activation(sigma)(x)
                x = layers.Dropout(rate=dropout)(x)
                if residual:
                    x = layers.add([x, y])
            return x

        return _feedforward, hid_dim[-1]
