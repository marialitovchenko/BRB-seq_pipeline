# Generated by Django 3.0.5 on 2020-05-19 13:26

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('userHome', '0003_auto_20200519_1010'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='project',
            name='sequencing',
        ),
        migrations.AddField(
            model_name='seqlibrary',
            name='project',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='userHome.Project'),
        ),
    ]
