from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='user-home'),
    path('tutorial/', views.tutorial, name='user-tutorial')
]
