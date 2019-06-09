
from fortranformat.FortranRecordReader import FortranRecordReader
from fortranformat.FortranRecordWriter import FortranRecordWriter

EL_FORMAT = '(I2,2F10.5,A4,I2)'
EL_READER = FortranRecordReader(EL_FORMAT)
EL_WRITER = FortranRecordWriter(EL_FORMAT)
IN1_UNIT_READER_v171 = FortranRecordReader('(20X,I1,2F10.1,I6)')
IN1_UNIT_READER_v142 = FortranRecordReader('(20X,I1,F7.1,F10.1,I6)')
VSP_VNS_READER = FortranRecordReader('(3X, 4ES19.12)')
STRUCT_LATT_PARAM_READER = FortranRecordReader('(6F10.6)')
STRUCT_LATT_PARAM_WRITER = FortranRecordWriter('(6F10.6)')

DEFAULT_RMTS = {
    "O": 1.4, 
    "Ti": 2.0,
}
DEFAULT_RMT = 1.8

DEFAULT_R0S = {
}
DEFAULT_R0 = 1e-4

DEFAULT_NPT = 781
