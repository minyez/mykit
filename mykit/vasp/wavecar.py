# conding = utf-8

class wavecarError(Exception):
    pass


class Wavecar:

    def __init__(self, WavFile='WAVECAR', verbose=True):
        
        self.filename = WavFile
        self.verbose = verbose
        with open(WavFile,'rb') as h_WavFile:
            #self.wavlines = h_WavFile.readlines()
            pass
