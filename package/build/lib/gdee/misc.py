"""
"""
import re
import json
import numpy as np


def _jsonfy(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()

    raise TypeError()


class DataContainer(dict):
    def __setattr__(self, name, value):
        self[name] = value

    def __getattr__(self, name):
        if name in self:
            return self[name]
        raise AttributeError(name)

    def jsonfy(self):
        return json.dumps(self, default=_jsonfy)


def get_valid_filename(name):
    name = str(name).strip().replace(" ", "_")
    return re.sub(r"(?u)[^-\w.]", "", name)
