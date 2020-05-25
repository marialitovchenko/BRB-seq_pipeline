from django import forms
from .models import Project, SeqLibrary, Specie, GenomeVersion, TrimGaloreParams

class SeqLibraryForm(forms.ModelForm):
    class Meta:
        model = SeqLibrary
        fields = ('RunID', 'LibraryID', 'SampleID', 'specie', 'genome')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['genome'].queryset = GenomeVersion.objects.none()
        if 'specie' in self.data:
            try:
                specie_id = int(self.data.get('specie'))
                self.fields['genome'].queryset = GenomeVersion.objects.filter(specie_id=specie_id).order_by('version')
            except (ValueError, TypeError):
                pass  # invalid input from the client; ignore and fallback to empty City queryset
        elif self.instance.pk:
            self.fields['city'].queryset = self.instance.specie.city_set.order_by('version')

class TrimGaloreParamsForm(forms.ModelForm):
    class Meta(object):
            model = TrimGaloreParams
            fields = '__all__'
            exclude = ('project',)
