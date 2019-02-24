# coding = utf-8
'''
'''
class programError(Exception):
    pass

class program:

    __progName = "n a"

    def __init__(self, program, **kwargs):
        self.__progName = program.lower()

    @property
    def progName(self):
        return self.__progName