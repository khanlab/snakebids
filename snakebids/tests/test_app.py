from pytest_mock.plugin import MockerFixture
from ..app import SnakeBidsApp
from .mock.config import config
from unittest.mock import patch
import pytest
import sys, copy

def init_snakebids_app(self):
    self.configfile_path="mock/config.yaml"
    self.snakefile = "mock/Snakefile"
    self.config = copy.deepcopy(config)



@pytest.fixture
def app(mocker: MockerFixture):
    mocker.patch.object(SnakeBidsApp, '__init__', init_snakebids_app)
    return SnakeBidsApp()


class TestCreateParser:
    mock_args_special= ["--derivatives", "path/to/nowhere"]
    mock_basic_args = ["script_name", "path/to/input", "path/to/output", "participant"]
    mock_all_args = mock_basic_args + mock_args_special 
    def test_snakebids_app_is_properly_mocked(self, app):
        assert isinstance(app, SnakeBidsApp)
        assert not hasattr(app, "parser")

    def test_fails_if_missing_arguments(self, app, mocker: MockerFixture):
        mocker.patch.object(sys, 'argv', self.mock_args_special)
        parser = app._SnakeBidsApp__create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args()

    def test_succeeds_if_given_positional_args(self, app, mocker: MockerFixture):
        mocker.patch.object(sys, 'argv', self.mock_basic_args)
        parser = app._SnakeBidsApp__create_parser()
        print(sys.argv)
        print(parser.parse_args())
        assert False
    
