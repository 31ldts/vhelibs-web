# -*- coding: utf-8 -*-
import logging
import os

from flask import Flask


def create_app(cache_dir=None):
    app = Flask(__name__, template_folder="templates", static_folder="static")
    app.config["SECRET_KEY"] = os.urandom(24)
    app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0  # never cache static files (dev convenience)
    app.config["CACHE_DIR"] = cache_dir or os.path.join(
        os.path.expanduser("~"), ".cache", "vhelibs"
    )
    os.makedirs(app.config["CACHE_DIR"], exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    with app.app_context():
        from app import routes  # noqa: F401
        app.register_blueprint(routes.bp)

    return app
