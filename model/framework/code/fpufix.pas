unit fpufix;
{ Mask floating-point exceptions so checkmol's geometry math does not trap
  with EInvalidOp on Linux (FPC unmasks FP exceptions by default). The
  initialization section runs before the checkmol program body. }
interface
implementation
uses Math;
initialization
  SetExceptionMask([exInvalidOp, exDenormalized, exZeroDivide,
                    exOverflow, exUnderflow, exPrecision]);
end.
