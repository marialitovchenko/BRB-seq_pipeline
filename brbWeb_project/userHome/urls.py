from django.urls import path
from .views import (ProjectListView, 
	ProjectDetailView, 
	ProjectCreateView,
	ProjectUpdateView,
	ProjectDeleteView,
    SeqLibraryCreateView)
from . import views

urlpatterns = [
    path('userHome/', ProjectListView.as_view(), name = 'user-home'),
    path('tutorial/', views.tutorial, name = 'user-tutorial'),
    path('project/<int:pk>/', ProjectDetailView.as_view(), 
        name = 'project-detail'),
    path('project/<int:pk>/update/', ProjectUpdateView.as_view(), 
        name = 'project-update'),
    path('project/<int:pk>/delete/', ProjectDeleteView.as_view(), 
        name = 'project-delete'),
    path('project/new/', ProjectCreateView.as_view(), 
        name = 'project-create'),
    path('project/<int:pk>/newSeqLibrary', SeqLibraryCreateView.as_view(), 
        name = 'SeqLibrary-create'),

    path('ajax/load-genomes/', views.load_genomes, name='ajax_load_genomes'),
]
