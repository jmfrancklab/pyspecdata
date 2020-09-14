{% for item in methods %}
{%- if item not in inherited_members %}
    ~{{ name }}.{{ item }}
{%- endif %}
{%- endfor %}
