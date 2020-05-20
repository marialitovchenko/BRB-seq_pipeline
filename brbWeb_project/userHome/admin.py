from django.contrib import admin
from .models import Project, SeqLibrary, Specie, GenomeVersion

admin.site.register(Project)
admin.site.register(Specie)
admin.site.register(GenomeVersion)
admin.site.register(SeqLibrary)