class ConfigError(Exception):
    """Exception raised for errors with the Snakebids config."""

    def __init__(self, msg):
        self.msg = msg
        super().__init__()