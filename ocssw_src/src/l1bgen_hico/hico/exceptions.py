class PVQException(ValueError):
    pass


class USGNCCoarseTimeException(PVQException):
    pass


class PVQInterpolation(PVQException):
    pass


class EmptyCSVException(PVQException):
    pass

class BadFormatCSVException(PVQException):
    pass
