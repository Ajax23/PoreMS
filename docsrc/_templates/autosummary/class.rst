.. raw:: html

    </div>
    <div class=col-md-9 content>

{{ objname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block members %}

   {% if members %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in members %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
