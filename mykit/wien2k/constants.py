
from fortranformat.FortranRecordReader import FortranRecordReader
from fortranformat.FortranRecordWriter import FortranRecordWriter

EL_FORMAT = '(I2,2F10.5,A4,I2)'
EL_READER = FortranRecordReader(EL_FORMAT)
EL_WRITER = FortranRecordWriter(EL_FORMAT)
IN1_UNIT_READER_v171 = FortranRecordReader('(20X,I1,2F10.1,I6)')
IN1_UNIT_READER_v142 = FortranRecordReader('(20X,I1,F7.1,F10.1,I6)')
