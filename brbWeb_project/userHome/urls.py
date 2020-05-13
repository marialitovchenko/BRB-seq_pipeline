from django.urls import path
from .views import (ProjectListView, 
	ProjectDetailView, 
	ProjectCreateView)
from . import views

urlpatterns = [
    path('userHome/', ProjectListView.as_view(), name = 'user-home'),
    path('project/<int:pk>/', ProjectDetailView.as_view(), name = 'project-detail'),
    path('project/new/', ProjectCreateView.as_view(), name = 'project-create'),
    path('tutorial/', views.tutorial, name = 'user-tutorial')
]
