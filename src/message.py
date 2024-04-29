import io
import logging
import os

class InputError(Exception):
    pass


logging.basicConfig(stream=io.StringIO(), level=logging.DEBUG)

# Create logger
logger = logging.getLogger()  # Main logger instance
logger.propagate=False

def configure(output_dir):
    if os.path.exists(os.path.join(output_dir,'warning.log')):
        os.remove(os.path.join(output_dir,'warning.log'))

    # Create handlers
    error_handler = logging.FileHandler(os.path.join(output_dir,'error.log'))
    warning_handler = logging.FileHandler(os.path.join(output_dir,'warning.log'))
    debug_handler = logging.FileHandler(os.path.join(output_dir,'debug.log'))
    
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



