# coding = utf8

from mykit.core.log import verbose

class controlError(Exception):
    pass

class control_map:

    @staticmethod
    def _tags_mapping(mapDict, progFrom, progTo, *tags, getAll=False):
        '''
        Args:
            m (dict): dictionary, with key-value pair of `'n a' tag name: 'n a' map (dict)`
            progFrom (str)
            progTo (str)
        '''
        if len(tags) == 0 and not getAll:
            return tuple()
        assert isinstance(mapDict, dict)
        for _v in mapDict.values():
            assert isinstance(_v, dict)
        verbose.print_cm_log("Use mapDict:", mapDict, level=3, depth=1)
        verbose.print_cm_log("To map tags:", tags, level=3, depth=1)
        _pF = progFrom
        _pT = progTo
        _d = {}
        for _map in mapDict.values():
            _d.update({_map.get(_pF, None): _map.get(_pT, None)})
        if None in _d:
            _d.update({None:None})
        if getAll:
            _d.pop(None, None)
            return tuple(_d.values())
        if len(tags) == 1:
            return (_d.get(tags[0], None),)
        return tuple(_d.get(_t, None) for _t in tags)

    @staticmethod
    def _val_mapping(m, progFrom, progTo, *val, getAll=False):
        raise NotImplementedError