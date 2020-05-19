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

class SeqLibrary(models.Model):
	project = models.ForeignKey(Project, related_name = 'seqData', null = True,
		on_delete = models.CASCADE)
	RunID = models.CharField(max_length = 100)
	LibraryID = models.CharField(max_length = 100)
	SampleID = models.CharField(max_length = 100)
	specieChoices = [ ('H.sapiens', 'H.sapiens'),
	('M.musculus', 'M.musculus'),
	('D.melanogaster', 'D.melanogaster')]
	Specie = models.CharField(max_length = 20, choices = specieChoices, 
		default = 'H.sapiens',)
	genomeChoices = [
    ('H.sapiens', (
            ('vinyl', 'Vinyl'),
            ('cd', 'CD'),
        )
    ),
    ('M.musculus', (
            ('vhs', 'VHS Tape'),
            ('dvd', 'DVD'),
        )
    ),
    ('unknown', 'Unknown'),
	]

	def __str__ (self):
		return self.RunID