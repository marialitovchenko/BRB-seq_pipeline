from django.shortcuts import render

def index(request):
    return render(request, 'public_website/index.html', {'title': 'Alithea Genomics'})