from pygments import token
from pygments.styles.dracula import DraculaStyle


class CustomDraculaStyle(DraculaStyle):

    styles = DraculaStyle.styles | {
        token.Generic.Output: "#f8f8f2",
    }
