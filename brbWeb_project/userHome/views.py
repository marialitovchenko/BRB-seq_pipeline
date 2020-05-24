import csv, io
from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.views.generic import (ListView, 
	DetailView, 
	CreateView,
	UpdateView,
	DeleteView)
from django_datatables_view.base_datatable_view import BaseDatatableView
from django.utils.html import escape
from .models import Project, SeqLibrary, Specie, GenomeVersion
from .forms import SeqLibraryForm

def home(request):
	context = {
		'projects' : Project.objects.all()
	}
	return render(request, 'userHome/home.html', context)

class ProjectListView(ListView):
	model = Project
	template_name = 'userHome/home.html'
	context_object_name = 'projects'
	ordering = ['-date_created']
	#paginate_by = 5

class ProjectDetailView(DetailView):
	model = Project

class ProjectCreateView(LoginRequiredMixin, CreateView):
	model = Project
	fields = ['name', 'description']

	def form_valid(self, form):
		form.instance.author = self.request.user
		render(request, 'SeqLibrary-create', context)
		#return super().form_valid(form)

class ProjectUpdateView(LoginRequiredMixin, UserPassesTestMixin, UpdateView):
	model = Project
	fields = ['name', 'description']

	def form_valid(self, form):
		form.instance.author = self.request.user
		return super().form_valid(form)

	def test_func(self):
		project = self.get_object()
		if self.request.user == project.author :
			return True
		return False

class ProjectDeleteView(LoginRequiredMixin, UserPassesTestMixin, DeleteView):
	model = Project
	success_url = '/userHome/userHome'

	def test_func(self):
		project = self.get_object()
		if self.request.user == project.author :
			return True
		return False

class SeqLibraryCreateView(LoginRequiredMixin, CreateView):
	model = SeqLibrary
	form_class = SeqLibraryForm

	def dispatch(self, request, *args, **kwargs):
		self.project = get_object_or_404(Project, pk=kwargs['project_pk'])
		return super().dispatch(request, *args, **kwargs)

	def form_valid(self, form):
		form.instance.project = self.project
		return super().form_valid(form)

def load_genomes(request):
    specie_id = request.GET.get('specie')
    genomes = GenomeVersion.objects.filter(specie_id=specie_id).order_by('version')
    return render(request, 'userHome/genomes_dropdown_form.html', {'genomes': genomes})

def SeqLibrary_upload(request, project_pk):
	template = "userHome/SeqLibrary_upload.html"
	promt = {
		'order' : 'Order of the CSV should be…'
	}

	if request.method == "GET":
		return render(request, template, promt)

	csv_file = request.FILES['file']

	if not csv_file.name.endswith('.csv'):
		messages.error(request, 'This is not a csv file')

	data_set = csv_file.read().decode('UTF-8')
	io_string = io.StringIO(data_set)
	next(io_string) # skip the header
	for column in csv.reader(io_string, delimiter=',', quotechar="|"):
		_,created = SeqLibrary.objects.update_or_create(
			project = get_object_or_404(Project, pk=project_pk),
			RunID = column[0],
			LibraryID = column[1],
			SampleID = column[2],
			specie = get_object_or_404(Specie, name = column[3]),
			genome = get_object_or_404(GenomeVersion, version = column[4])
		)
		
	context = {}
	return redirect('project-detail', project_pk)

def tutorial(request):
    return render(request, 'userHome/tutorial.html', {'title': 'Test title'})