from django.urls import path
from .views import (ProjectListView, 
	ProjectDetailView, 
	ProjectCreateView,
	ProjectUpdateView)
from . import views

urlpatterns = [
    path('userHome/', ProjectListView.as_view(), name = 'user-home'),
    path('project/<int:pk>/', ProjectDetailView.as_view(), name = 'project-detail'),
    path('project/<int:pk>/update/', ProjectUpdateView.as_view(), name = 'project-update'),
    path('project/new/', ProjectCreateView.as_view(), name = 'project-create'),
    path('tutorial/', views.tutorial, name = 'user-tutorial')
]
