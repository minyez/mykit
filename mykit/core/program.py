# coding = utf-8
'''
'''
class programError(Exception):
    pass

class program:

    __progName = "n a"

    def __init__(self, **progargs):
        if "program" in progargs:
            self.__progName = progargs["program"].lower()
        else:
            raise programError("No program name is specified by keyword `program`.")

    @property
    def progName(self):
        return self.__progName