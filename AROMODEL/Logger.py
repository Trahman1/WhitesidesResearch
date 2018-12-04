# class of loggers takes in module name

import datetime

class Logger():

	def __init__(self, name):
		self.name = name
		self.log(self.name + " initialized!")

	def log(self, str_to_log):
		with open ("log.txt", "a") as log:
			str_to_log = self.name + ": " + str_to_log + " " + str(datetime.datetime.today()) + "\n\n"
			print str_to_log
			log.write(str_to_log)

