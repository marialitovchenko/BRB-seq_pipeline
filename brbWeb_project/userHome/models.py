from django.db import models
from django.utils import timezone
from django.contrib.auth.models import User
from django.urls import reverse

class Project(models.Model):
	author = models.ForeignKey(User, on_delete = models.CASCADE)
	name = models.CharField(max_length = 100)
	description = models.TextField()
	date_created = models.DateTimeField(auto_now_add = True)
	date_updated = models.DateTimeField(auto_now = True)

	def __str__ (self):
		return self.name

	def get_absolute_url(self):
		return reverse('project-detail', kwargs = {'pk': self.pk})

class Specie(models.Model):
	name = models.CharField(max_length = 100)

	def __str__(self):
		return self.name

class GenomeVersion(models.Model):
	specie = models.ForeignKey(Specie, max_length = 100, 
		on_delete = models.CASCADE)
	version = models.CharField(max_length = 100)

	def __str__ (self):
		return self.version

class SeqLibrary(models.Model):
	project = models.ForeignKey(Project, related_name = 'seqData', null = True,
		on_delete = models.CASCADE)
	RunID = models.CharField(max_length = 100)
	LibraryID = models.CharField(max_length = 100)
	SampleID = models.CharField(max_length = 100)
	specie = models.ForeignKey(Specie, null = True, on_delete = models.SET_NULL)
	genome = models.ForeignKey(GenomeVersion, null = True, on_delete = models.SET_NULL)

	def __str__ (self):
		return self.RunID

	def get_absolute_url(self):
		return reverse('project-detail', kwargs = {'pk': self.project.pk})