# coding = utf8

class controlError(Exception):
    pass

class control:

    @classmethod
    def _tag_mapping(cls, map, progFrom, progTo, *tag, getAll=False):
        raise controlError