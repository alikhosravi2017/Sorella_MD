#/usr/bin/env python3


class Car(object):
	def __init__(self, age):
		super(Car, self).__init__()
		self.age = age
		self.age_sq = age**2
		
	def driving(self,speed):
		self.speed = speed
		if self.age>18:
			print("Yeah he is driving")
		else:
			print("too you to drive")
		if self.speed<=100:
			print("he drives so fast")
		else:
			print("he drives like crazy")
			self.drunk = True

	def is_drunk(self):
		if self.drunk:
			print("Ali is drunk")



Ali = Car(25)
Ali.driving(120)
Ali.is_drunk()
