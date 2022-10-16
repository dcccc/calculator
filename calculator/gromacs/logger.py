import logging,sys

date_fmt = '%Y-%m-%d %I:%M:%S'
fmt = "[%(asctime)s][%(levelname)s] %(message)s "

class Logger():
    def __init__(self):
        self.logger = logging.getLogger()

        self.logger.setLevel(level=logging.INFO)

        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(level=logging.INFO)

        file_handler = logging.FileHandler(filename='workflow.log', encoding="utf-8", delay = False)
        file_handler.setLevel(level=logging.INFO)

        formatter = logging.Formatter(fmt, datefmt=date_fmt)
        file_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter )    

        self.logger.addHandler(file_handler)
        self.logger.addHandler(stream_handler)