"""An helper for sample dimensions. """


class SampleDimensions:
    """
    A sample is caracterised by three parameters width, thickness and length all in meters.
    """

    def __init__(self, width=1e-6, thickness=1e-6, length=1e-6):
        self.width = width
        self.thickness = thickness
        self.length = length

    def __repr__(self):
        return f"{self.width} x {self.thickness} x {self.length}"

    def geometric_factor(self):
        """
        Returns the geometric factor, or alpha, of the sample.
        alpha = w * t / L
        """
        return self.width * self.thickness / self.length
