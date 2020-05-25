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

class TrimGaloreParams(models.Model):
	project = models.ForeignKey(Project, null = True, on_delete = models.CASCADE)
	phred33 = models.CharField(max_length = 100)
	phred64 = models.CharField(max_length = 100)
	fastqc = models.CharField(max_length = 100)
	fastqc_args = models.CharField(max_length = 100)
	adapter = models.CharField(max_length = 100)
	adapter2 = models.CharField(max_length = 100)
	illumina = models.CharField(max_length = 100)
	nextera = models.CharField(max_length = 100)
	small_rna = models.CharField(max_length = 100)
	consider_already_trimmed = models.CharField(max_length = 100)
	max_length = models.CharField(max_length = 100)
	stringency = models.CharField(max_length = 100)
	e = models.CharField(max_length = 100)
	gzip = models.CharField(max_length = 100)
	dont_gzip = models.CharField(max_length = 100)
	length = models.CharField(max_length = 100)
	max_n = models.CharField(max_length = 100)
	trim_n = models.CharField(max_length = 100)
	output_dir = models.CharField(max_length = 100)
	no_report_file = models.CharField(max_length = 100)
	suppress_warn = models.CharField(max_length = 100)
	clip_R1 = models.CharField(max_length = 100)
	clip_R2 = models.CharField(max_length = 100)
	three_prime_clip_R1 = models.CharField(max_length = 100)
	three_prime_clip_R2 = models.CharField(max_length = 100)
	nextseq = models.CharField(max_length = 100)
	path_to_cutadapt = models.CharField(max_length = 100)
	basename = models.CharField(max_length = 100)
	cores = models.CharField(max_length = 100)
	hardtrim5 = models.CharField(max_length = 100)
	hardtrim3 = models.CharField(max_length = 100)
	clock = models.CharField(max_length = 100)
	polyA = models.CharField(max_length = 100)
	rrbs = models.CharField(max_length = 100)
	non_directional = models.CharField(max_length = 100)
	keep = models.CharField(max_length = 100)
	paired = models.CharField(max_length = 100)
	trim1 = models.CharField(max_length = 100)
	retain_unpaired = models.CharField(max_length = 100)
	length_1 = models.CharField(max_length = 100)
	length_2 = models.CharField(max_length = 100)


class STARParams(models.Model):
	project = models.ForeignKey(Project, null = True, on_delete = models.CASCADE)
	runThreadN = models.CharField(max_length = 100)