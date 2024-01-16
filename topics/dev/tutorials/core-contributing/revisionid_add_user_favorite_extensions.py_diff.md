
> <solution-title>``...alembic/versions_gxy/2ad8047d652e_add_user_favorite_extentions_table.py``</solution-title>
> 
> Possible changes to revision template:
> 
> ```diff
> new file mode 100644
> index 0000000000..e8e5fe0ea3
> --- /dev/null
> +++ b/lib/galaxy/model/migrations/alembic/versions_gxy/2ad8047d652e_add_user_favorite_extentions_table.py
> @@ -0,0 +1,29 @@
> +"""Add user_favorite_extentions table
> +
> +Revision ID: 2ad8047d652e
> +Revises: 186d4835587b
> +Create Date: 2022-07-07 00:44:21.992162
> +
> +"""
> +from alembic import op
> +import sqlalchemy as sa
> +
> +
> +# revision identifiers, used by Alembic.
> +revision = '2ad8047d652e'
> +down_revision = '186d4835587b'
> +branch_labels = None
> +depends_on = None
> +
> +
> +def upgrade():
> +    op.create_table(
> +        'user_favorite_extension',
> +        sa.Column('id', sa.Integer, primary_key=True),
> +        sa.Column('user_id', sa.Integer, sa.ForeignKey("galaxy_user.id")),
> +        sa.Column('value', sa.String),
> +    )
> +
> +
> +def downgrade():
> +    op.drop_table('user_favorite_extension')
> ```
{: .solution }
