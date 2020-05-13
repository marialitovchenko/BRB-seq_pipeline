from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import (ListView, 
	DetailView, 
	CreateView)
from .models import Project

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

class ProjectDetailView(DetailView):
	model = Project

class ProjectCreateView(CreateView):
	model = Project
	fields = ['name', 'description']

def tutorial(request):
    return render(request, 'userHome/tutorial.html', {'title': 'Test title'})