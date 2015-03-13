
#ifndef MOUSE_H
#define MOUSE_H

#include<QVector2D>
#include<QPoint>

class Mouse {
	private:
		QPoint oldPosition;
		QPoint newPosition;
		QPoint _pressedCell;
		QVector2D _deltaPosition;

	public:
		void setPosition(const QPoint& _pos);
		QPoint& pressedCell(const QVector2D& _cellPixelSize);
		QVector2D& deltaPosition();
};

#endif