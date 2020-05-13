from django.urls import path
from .views import ProjectListView, ProjectDetailView
from . import views

urlpatterns = [
    path('userHome/', ProjectListView.as_view(), name = 'user-home'),
    path('project/<int:pk>/', ProjectDetailView.as_view(), name = 'project-detail'),
    path('tutorial/', views.tutorial, name = 'user-tutorial')
]
