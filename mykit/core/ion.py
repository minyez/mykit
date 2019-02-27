# coding = utf-8

from mykit.core._control import tags_mapping, prog_mapper
from mykit.core.log import verbose

class ionError(Exception):
    pass


class ion_control(verbose):

    _tagMaps = {
                "nSteps": {"mykit": "nSteps"}
               }
    _ionTagMaps = _tagMaps

    def __init__(self, **iontags):
        pass

    # @classmethod
    # def map_tags(cls, *tags, progFrom="mykit", progTo="mykit", getAll=False):
    #     _mapDict = {"sometag": {"mykit": "sometag"}}
    #     return tags_mapping(_mapDict, progFrom, progTo, *tags, getAll=getAll)