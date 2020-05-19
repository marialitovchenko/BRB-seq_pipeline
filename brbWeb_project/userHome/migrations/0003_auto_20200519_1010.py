# Generated by Django 3.0.5 on 2020-05-19 10:10

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('userHome', '0002_seqlibrary'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='seqlibrary',
            name='project',
        ),
        migrations.AddField(
            model_name='project',
            name='sequencing',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='userHome.SeqLibrary'),
        ),
    ]