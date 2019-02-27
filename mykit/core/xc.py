# coding = utf8
'''
'''
from mykit.core.log import verbose
from mykit.core._control import tags_mapping, prog_mapper

class XCError(Exception):
    pass


class xc_control(verbose, prog_mapper):
    '''base class that controls exchange-correlation setup
    '''

    _tagMaps = {
                "gga" : {"mykit":"gga", "vasp":"GGA"},
                "metagga" : {"mykit":"metagga", "vasp":"METAGGA"},
                "percentGgaCorr": {"mykit": "percentGgaCorr", "vasp":"AGGAC"},
                "percentGgaExch": {"mykit": "percentGgaExch", "vasp":"AGGAX"},
                "percentExactExch": {"mykit": "percentExactExch", "vasp":"AExX"},
                "percentLdaCorr": {"mykit": "percentLdaCorr", "vasp":"ALDAC"},
               }
    _xcTagMaps = _tagMaps
    _xcValMaps = {}
    _xcTags = {}

    def __init__(self, progName, **xctags):
        self._parse_xctags(progName, **xctags)

    def parse_tags(self, progName, **xctags):
        self._parse_xctags(progName, **xctags)

    def _parse_xctags(self, progName, **xctags):
        if len(xctags) == 0:
            return
        self.print_log(" In _parse_xctags. Parsing: ", xctags, depth=1, level=3)
        for _origTag, _v in xctags.items():
            if _origTag == None:
                continue
            elif _origTag in self._xcTagMaps:
                self._xcTags.update({_origTag:_v})
            else:
                for _xct, _xcmap in self._xcTagMaps.items():
                    if _origTag == _xcmap.get(progName, None):
                        self._xcTags.update({_xct: _v})
                        break
        self.print_log("End _parse_xctags, now xcTags: ", self._xcTags, depth=1, level=3)
        
    def delete_tags(self, progName, *tags):
        self._pop_xctags(progName, *tags)

    def pop_tags(self, progName, *tags):
        return self._pop_xctags(progName, *tags)

    def _pop_xctags(self, progName, *tags):
        _vals = self._xctag_vals(progName, *tags, delete=True)
        return _vals

    def _get_one_xctag(self, xcTagName):
        return self._xcTags.get(xcTagName, None)

    def tag_vals(self, progName, *tags):
        '''find values of xc tags from program-specific tags of progName
        
        The tags name depends on the program, i.e. ``progName``.

        Note:
            Tags that belong to xcTagMaps will also get their value instead of None,
            even when progName is not assigned, i.e. use "mykit"

        Args:
            tags (str): names of tags to request values
            progName (str): the program of which the tags belong to.

        Returns:
            list containing all tag values.
        '''
        return self._xctag_vals(progName, *tags)

    def _xctag_vals(self, progName, *tags, delete=False):
        if len(tags) == 0:
            return []
        self.print_log("In _xctag_vals, search {} tags of {}: ".format(len(tags), progName), tags, level=3, depth=1)
        _xctags = xc_control.map_to_mykit_tags(*tags, progFrom=progName)
        # self.print_log("_xctags:", _xctags, level=3, depth=2)
        _vals = list(map(self._get_one_xctag, _xctags))
        # self.print_log("_values:", _vals, level=3, depth=2)
        if progName != "mykit":
            for _i, _v in enumerate(_vals):
                if _v == None:
                    if tags[_i] in self._xcTags:
                        _vals[_i] = self._xcTags[tags[_i]]
        if delete:
            for _i, _v in enumerate(_xctags):
                if _v in self._xcTags:
                    del self._xcTags[_v]
                if tags[_i] in self._xcTags:
                    del self._xcTags[tags[_i]]
        self.print_log("Found values:", _vals, level=3, depth=3)
        return _vals

    @classmethod
    def map_tags(cls, *tags, progFrom="mykit", progTo="mykit", getAll=False):
        _pF = progFrom.lower()
        _pT = progTo.lower()
        return tags_mapping(cls._xcTagMaps, _pF, _pT, *tags, getAll=getAll)

    @property
    def xcTags(self):
        return self._xcTags