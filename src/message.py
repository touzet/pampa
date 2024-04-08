import io
import logging


class InputError(Exception):
    pass

logging.basicConfig(stream=io.StringIO(), level=logging.DEBUG)

# Create loggers
logger = logging.getLogger()  # Main logger instance
logger.propagate=False

# Create handlers
error_handler = logging.FileHandler('error.log')
warning_handler = logging.FileHandler('warning.log')
debug_handler = logging.FileHandler('debug.log')

# Set levels for handlers
error_handler.setLevel(logging.ERROR)  
warning_handler.setLevel(logging.WARNING)
debug_handler.setLevel(logging.DEBUG)

#Create formatters and set them for handlers
formatter = logging.Formatter('')
error_handler.setFormatter(formatter)
warning_handler.setFormatter(formatter)
debug_handler.setFormatter(formatter)

# Add handlers to the logger
logger.addHandler(error_handler)
logger.addHandler(warning_handler)
logger.addHandler(debug_handler)

def escape(s):
    logger.error("\n"+s+"\n")
    raise InputError

def warning(s):
    logger.warning(s)

def debug(s):
    logger.debug(s)



