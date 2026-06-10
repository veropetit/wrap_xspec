import xspec
# Configure XSpec globally for the entire package session
xspec.Xset.allowPrompting = False

from . import plot
from .model import *
from .SpectrumLoader import SpectrumLoader

# This tells Pylance: "These are public. Expose their full types and docstrings."
__all__ = [
    "SpectrumLoader",
    "plot"
]