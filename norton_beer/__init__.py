try:
    from nbapo.version import COMMIT as __commit__
    from nbapo.version import REVISION as __revision__
except ImportError:
    __version__ = "unbuilt-dev"
    __revision__ = "unbuilt-dev"
    __commit__ = "unbuilt-dev"


def get_info():
    return (f"REV{__revision__}_COM{__commit__}")
