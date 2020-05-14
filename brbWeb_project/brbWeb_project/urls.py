from django.contrib import admin
from django.contrib.auth import views as auth_views
from django.urls import path, include
from django.conf import settings
from django.conf.urls.static import static
from users import views as user_view

urlpatterns = [
    path('admin/', admin.site.urls),
    path('register/', user_view.register, name = 'Register'),
    path('login/', 
         auth_views.LoginView.as_view(template_name = 'users/login.html'), 
         name = 'Log in'),
    path('logout/',
         auth_views.LogoutView.as_view(template_name = 'users/logout.html'), 
         name = 'Log out'),
    path('password-reset/',
         auth_views.PasswordResetView.as_view(template_name = 'users/password_reset.html'), 
         name = 'password reset'),
    path('password-reset/done/',
         auth_views.PasswordResetDoneView.as_view(template_name = 'users/password_reset_done.html'), 
         name = 'password reset done'),
    path('password-reset-confirm/<uidb64>/<token>/',
         auth_views.PasswordResetConfirmView.as_view(template_name = 'users/password_reset_confirm.html'), 
         name = 'password_reset_confirm'),
    path('profile/', user_view.profile, name = 'Profile'),
    path('userHome/', include('userHome.urls')),
    path('', include('public_website.urls')),
]

if settings.DEBUG: 
    urlpatterns += static(settings.MEDIA_URL, document_root = settings.MEDIA_ROOT)