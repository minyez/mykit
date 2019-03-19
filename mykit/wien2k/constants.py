
from fortranformat.FortranRecordReader import FortranRecordReader
from fortranformat.FortranRecordWriter import FortranRecordWriter

EXCEPTION_FORMAT = '(I2,2F10.5,A4,I2)'
EXCEPTION_READER = FortranRecordReader(EXCEPTION_FORMAT)
EXCEPTION_WRITER = FortranRecordWriter(EXCEPTION_FORMAT)
IN1_UNIT_READER_v171 = FortranRecordReader('(20X,I1,2F10.1,I6)')
IN1_UNIT_READER_v142 = FortranRecordReader('(20X,I1,F7.1,F10.1,I6)')
