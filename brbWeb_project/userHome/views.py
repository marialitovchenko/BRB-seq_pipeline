from django.shortcuts import render
from django.http import HttpResponse
from .models import Project

def home(request):
	context = {
		'projects' : Project.objects.all()
	}
	return render(request, 'userHome/home.html', context)

def tutorial(request):
    return render(request, 'userHome/tutorial.html', {'title': 'Test title'})