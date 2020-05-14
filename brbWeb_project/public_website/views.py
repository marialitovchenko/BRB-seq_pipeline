from django.shortcuts import render

def public_home(request):
    return render(request, 'public_website/public_home.html', {'title': 'Alithea Genomics'})