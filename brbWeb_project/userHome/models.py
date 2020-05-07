from django.db import models
from django.utils import timezone
from django.contrib.auth.models import User

class Project(models.Model):
	author = models.ForeignKey(User, on_delete = models.CASCADE)
	name = models.CharField(max_length = 100)
	description = models.TextField()
	date_created = models.DateTimeField(auto_now_add = True)
	date_updated = models.DateTimeField(auto_now = True)

	def __str__ (self):
		return self.name