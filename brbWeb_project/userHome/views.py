from django.shortcuts import render
from django.http import HttpResponse
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.views.generic import (ListView, 
	DetailView, 
	CreateView,
	UpdateView,
	DeleteView)
from .models import Project, SeqLibrary

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
	paginate_by = 5

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
	fields = ['RunID', 'LibraryID', 'SampleID', 'Specie']

	def form_valid(self, form):
		form.instance.author = self.request.project
		return super().form_valid(form)

def tutorial(request):
    return render(request, 'userHome/tutorial.html', {'title': 'Test title'})